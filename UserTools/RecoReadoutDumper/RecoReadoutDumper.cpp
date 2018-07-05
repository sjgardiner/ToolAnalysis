// ToolAnalysis includes
#include "DataModel.h"
#include "RecoReadoutDumper.h"

// standard library includes
#include <map>
#include <memory>
#include <string>
#include <vector>

RecoReadoutDumper::RecoReadoutDumper() : Tool()
{}

bool RecoReadoutDumper::Initialise(const std::string config_file,
  DataModel& data)
{
  // Load configuration file variables
  if ( !config_file.empty() ) m_variables.Initialise(config_file);

  // Assign transient data pointer
  m_data = &data;

  // Prepare to write data to the output ROOT file
  TTree* out_tree_ = nullptr;

  std::string output_file_name;
  m_variables.Get("OutputFile", output_file_name);

  out_file_ = std::unique_ptr<TFile>( new TFile(output_file_name.c_str(),
    "recreate") );

  annie::RecoReadout* dummy_rr = nullptr;

  out_tree_ = new TTree("reco_readouts", "dumped RecoReadout objects");
  out_tree_->Branch("reco_readout", "annie::RecoReadout", &dummy_rr);
  //out_tree_->Branch("hefty_trigger_mask", &hefty_trigger_mask_,
  //  "hefty_trigger_mask/I");
  out_tree_->Branch("ncv_position", &ncv_position_, "ncv_position/I");
  out_tree_->Branch("hefty_mode", &hefty_mode_, "hefty_mode/O");
  out_tree_->Branch("run_number", &run_number_, "run_number/I");
  out_tree_->Branch("subrun_number", &subrun_number_, "subrun_number/I");
  out_tree_->Branch("event_number", &event_number_, "event_number/I");

  return true;
}

bool RecoReadoutDumper::Execute() {

  int verbosity;
  m_variables.Get("verbose", verbosity);

  // Get a pointer to the ANNIEEvent Store
  auto* annie_event = m_data->Stores["ANNIEEvent"];

  if (!annie_event) {
    Log("Error: The RecoReadoutDumper tool could not find the ANNIEEvent Store",
      0, verbosity);
    return false;
  }

  // Load the labels describing the type of data stored in each minibuffer
  std::vector<MinibufferLabel> mb_labels;

  get_object_from_store("MinibufferLabels", mb_labels, *annie_event);
  check_that_not_empty("MinibufferLabels", mb_labels);

  // Decide whether we're using Hefty vs. non-Hefty data by checking whether
  // the HeftyInfo object is present
  hefty_mode_ = annie_event->Has("HeftyInfo");

  // One of these objects will be used to get minibuffer timestamps
  // depending on whether we're using Hefty mode or non-Hefty mode.
  HeftyInfo hefty_info; // used for Hefty mode only
  std::vector<TimeClass> mb_timestamps; // used for non-Hefty mode only
  size_t num_minibuffers = 0u;

  if (hefty_mode_) {
    get_object_from_store("HeftyInfo", hefty_info, *annie_event);
    num_minibuffers = hefty_info.num_minibuffers();

    if ( num_minibuffers == 0u ) {
      Log("Error: The RecoReadoutDumper tool found an empty HeftyInfo entry", 0,
        verbosity);
      return false;
    }

    // Exclude beam spills (or source triggers) near the end of a full
    // multi-minibuffer readout that included extra self-triggers in the Hefty
    // window that could not be recorded. This is indicated in the heftydb
    // TTree via More[39] == 1 and in the HeftyInfo object by more() == true.
    if ( hefty_info.more() ) {
      // Find the first beam or source minibuffer counting backward from the
      // end of the full readout
      size_t mb = num_minibuffers - 1u;
      for (; mb > 0u; --mb) {
        MinibufferLabel label = mb_labels.at(mb);
        if ( label == MinibufferLabel::Beam
          || label == MinibufferLabel::Source ) break;
      }
      // Exclude the minibuffers from the incomplete beam or source trigger's
      // Hefty window by setting a new value of num_minibuffers. This will
      // prematurely end the loop over minibuffers below.
      num_minibuffers = mb;
    }
  }
  else {
    // non-Hefty data
    get_object_from_store("MinibufferTimestamps", mb_timestamps, *annie_event);
    check_that_not_empty("MinibufferTimestamps", mb_timestamps);
    num_minibuffers = mb_timestamps.size();
    // Trigger masks are not saved in the tree for non-Hefty mode
    hefty_trigger_mask_ = 0;
  }

  // Load the beam status objects for each minibuffer
  std::vector<BeamStatus> beam_statuses;
  get_object_from_store("BeamStatuses", beam_statuses, *annie_event);
  check_that_not_empty("BeamStatuses", beam_statuses);

  // Load run, subrun, and event numbers
  get_object_from_store("RunNumber", run_number_, *annie_event);
  get_object_from_store("SubRunNumber", subrun_number_, *annie_event);
  get_object_from_store("EventNumber", event_number_, *annie_event);

  // Determine the NCV position based on the run number
  ncv_position_ = get_NCV_position(run_number_);

  // Load the reconstructed ADC hits
  std::map<ChannelKey, std::vector< std::vector<ADCPulse> > > adc_hits;

  get_object_from_store("RecoADCHits", adc_hits, *annie_event);
  check_that_not_empty("RecoADCHits", adc_hits);

  // Create an empty annie::RecoReadout object
  annie::RecoReadout* temp_rr = new annie::RecoReadout( event_number_ );

  for (const auto& pair : adc_hits) {

    const auto& channel_key = pair.first;
    const auto& minibuffer_pulses = pair.second;

    int pulse_pmt_id = channel_key.GetDetectorElementIndex();

    // Flag that vetos minibuffers because the last beam minibuffer failed
    // the quality cuts.
    bool beam_veto_active = false;

    for (size_t mb = 0; mb < num_minibuffers; ++mb) {

      MinibufferLabel event_mb_label = mb_labels.at(mb);

      // If this is Hefty mode data, save the trigger mask for the
      // current minibuffer. This is distinct from the MinibufferLabel
      // assigned to the event_label_ variable above. The trigger
      // mask branch isn't used for non-Hefty data.
      if ( hefty_mode_ ) hefty_trigger_mask_ = hefty_info.label(mb);
      else hefty_trigger_mask_ = 0;

      // BEAM QUALITY CUT
      // Skip beam minibuffers with bad or missing beam status information
      // TODO: consider printing a warning message here
      const auto& beam_status = beam_statuses.at(mb);
      const auto& beam_condition = beam_status.condition();
      if (beam_condition == BeamCondition::Missing
        || beam_condition == BeamCondition::Bad)
      {
        // Skip all beam and Hefty window minibuffers until a good-quality beam
        // spill is found again
        beam_veto_active = true;
      }
      if ( beam_veto_active && beam_condition == BeamCondition::Ok ) {
        // We've found a new beam minibuffer that passed the quality check,
        // so disable the beam quality veto
        beam_veto_active = false;
      }
      if (beam_veto_active && (event_mb_label == MinibufferLabel::Hefty
          || event_mb_label == MinibufferLabel::Beam))
      {
        // Bad beam minibuffers and Hefty window minibuffers belonging to the
        // bad beam spill need to be skipped. Since other minibuffers (e.g.,
        // cosmic trigger minibuffers) are not part of the beam "macroevent,"
        // they may still be processed normally.
        continue;
      }

      for (const ADCPulse& pulse : minibuffer_pulses.at(mb) ) {
        // For non-Hefty mode, the neutron capture candidate event time
        // is simply its timestamp relative to the start of the single
        // minibuffer.
        pulse_start_time_ns_ = pulse.start_time().GetNs();

        // For Hefty mode, the pulse time within the current minibuffer
        // needs to be added to the TSinceBeam value for Hefty window
        // minibuffers (labeled by the RawLoader tool with
        // MinibufferLabel::Hefty). This should be zero for all other
        // minibuffers (it defaults to that value in the heftydb TTree).
        //
        // NOTE: event times for minibuffer labels other than "Beam", "Source",
        // and "Hefty" are calculated relative to the start of the minibuffer,
        // not the beam, source trigger, etc. Only make timing plots using
        // those 3 labels unless you're interested in single-minibuffer timing!
        if ( hefty_mode_ && event_mb_label == MinibufferLabel::Hefty) {
          // The name "TSinceBeam" is used in the heftydb tree for the time
          // since a beam *or* a source trigger, since the two won't be used
          // simultaneously.
          pulse_start_time_ns_ += hefty_info.t_since_beam(mb);
        }

        auto card_channel_pair = pmt_ID_and_card_channel_bimap
          .left.at(pulse_pmt_id);

        // TODO: add real peak time here (not just zero it out like we do
        // currently)
        temp_rr->add_pulse(card_channel_pair.first, card_channel_pair.second,
          mb, annie::RecoPulse(pulse_start_time_ns_, 0u, pulse.baseline(),
          pulse.sigma_baseline(), pulse.raw_area(), pulse.raw_amplitude(),
          pulse.amplitude(), pulse.charge()));


        Log("Found pulse on channel " + std::to_string(pulse_pmt_id)
          + " in run " + std::to_string(run_number_) + " subrun "
          + std::to_string(subrun_number_) + " event "
          + std::to_string(event_number_) + " in minibuffer "
          + std::to_string(mb) + " at "
          + std::to_string(pulse_start_time_ns_) + " ns", 3, verbosity);
      }
    }
  }

  out_tree_->SetBranchAddress("reco_readout", &temp_rr);
  out_tree_->Fill();

  return true;
}

int RecoReadoutDumper::get_NCV_position(uint32_t run_number) const {
  if (run_number >= 635u && run_number < 704u) return 1;
  if (run_number >= 704u && run_number < 802u) return 2;
  if (run_number >= 802u && run_number < 808u) return 3;
  if (run_number >= 808u && run_number < 813u) return 4;
  if (run_number == 813u) return 5;
  if (run_number == 814u) return 6;
  if (run_number >= 815u && run_number < 825u) return 7;
  if (run_number >= 825u && run_number < 883u) return 8;
  else return 0;
}

bool RecoReadoutDumper::Finalise() {
  out_tree_->Write();
  out_file_->Write();
  out_file_->Close();
  return true;
}
