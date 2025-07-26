#ifndef __Stntuple_gui_TStnWidgetID_hh__
#define __Stntuple_gui_TStnWidgetID_hh__
//-----------------------------------------------------------------------------

enum GuiCommandID {
  M_FILE_OPEN,
  M_FILE_SAVE,
  M_FILE_SAVEAS,
  M_FILE_EXIT,

  M_TEST_DLG,
  M_TEST_MSGBOX,
  M_TEST_SLIDER,
  M_TEST_SHUTTER,
  M_TEST_PROGRESS,

  M_EDIT_EDITOR,
  M_EDIT_UNDO,
  M_EDIT_CLEARPAD,
  M_EDIT_CLEARCANVAS,

  M_OPTION_EVENT_STATUS,
  M_OPTION_AUTO_EXEC,
  M_OPTION_AUTO_RESIZE,
  M_OPTION_RESIZE_CANVAS,
  M_OPTION_MOVE_OPAQUE,
  M_OPTION_RESIZE_OPAQUE,
  M_OPTION_REFRESH,
  M_OPTION_STATISTICS,
  M_OPTION_HIST_TITLE,
  M_OPTION_FIT_PARAMS,
  M_OPTION_CAN_EDIT,

  M_HELP_CONTENTS,
  M_HELP_SEARCH,
  M_HELP_ABOUT,

  M_OPEN_XY,
  M_OPEN_RZ,
  M_OPEN_TZ,
  M_OPEN_PHIZ,
  M_OPEN_CAL,
  M_OPEN_CRV,
  M_OPEN_VST,
  

  M_PRINT_STRAW_H,
  M_PRINT_COMBO_H,

  M_OPEN_LEGEND_XY, // legend popup
  M_OPEN_LEGEND_CAL,
  M_OPEN_LEGEND_CRV,
  M_PRINT_TO_PNG // printing functionality

};

enum GuiButtonID {
					// checkbutton ID's
  kDisplayHelices       =  5,
  kDisplayTracks        =  6,
  kDisplaySimParticles  =  7,
  kDisplayOnlyTCHits    =  8,
  kDisplaySH            =  9,
  kIgnoreComptonHits    = 10,
  kIgnoreProtonHits     = 11,
  kIgnoreProtons        = 12,
					// not sure what these are... but ...
                                        // they are needed for something
  M_DISPLAY_SH          = 21,
  M_DISPLAY_CH          = 22,
  M_IGNORE_COMPTON_HITS = 23,
  M_IGNORE_PROTON_HITS  = 24,
  M_IGNORE_PROTONS      = 25

};

#endif

