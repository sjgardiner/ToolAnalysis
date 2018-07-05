// RecoReadoutDumper tool
//
// Loads objects in the ANNIEEvent store with data from an event in an
// ANNIE Phase I raw data file. Uses code borrowed from the previous recoANNIE
// framework.
//
// Steven Gardiner <sjgardiner@ucdavis.edu>
#pragma once

// standard library includes
#include <iostream>
#include <memory>
#include <string>

// ToolAnalysis includes
#include "Tool.h"

// recoANNIE includes
#include "RecoReadout.h"
#include "RecoPulse.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

class RecoReadoutDumper : public Tool {

  public:

    RecoReadoutDumper();
    bool Initialise(const std::string config_file, DataModel& data) override;
    bool Execute() override;
    bool Finalise() override;

    int get_NCV_position(uint32_t run_number) const;

  protected:

    template <typename T, typename AStore> bool get_object_from_store(
      const std::string& object_label, T& obj, AStore& s)
    {
      Log("Retrieving \"" + object_label + "\" from a Store", 4, 0);
      bool got_object = s.Get(object_label, obj);

      // Check for problems
      if ( !got_object ) {
        Log("Error: The RecoReadoutDumper tool could not find the "
          + object_label + " entry", 0, 0);
        return false;
      }

      return true;
    }

    template <typename T> inline auto check_that_not_empty(
      const std::string& object_label, T& obj)
      -> decltype( std::declval<T&>().empty() )
    {
      bool is_empty = obj.empty();
      if ( is_empty ) {
        Log("Error: The RecoReadoutDumper tool found an empty " + object_label
          + " entry", 0, 0);
      }

      return !is_empty;
    }

    bool hefty_mode_ = false;
    int hefty_trigger_mask_ = 0;
    int ncv_position_ = 0;
    int run_number_ = 0;
    int subrun_number_ = 0;
    int event_number_ = 0;
    size_t pulse_start_time_ns_ = 0u;

    std::unique_ptr<TFile> out_file_;
    TTree* out_tree_ = nullptr;
};
