//-----------------------------------------------------------------------------
// available job configurations:
// ----------------------------
//  lum               :
//
// env variables used   :
// ----------------------
//----------------------------------------------------------------------------- 
// run range for the publication: 141544-156487 (before the winter'2003 shutdown)
//                                good run list: "WZ_PRD"
// last run before the Sept'2003 shutdown: 168899
//-----------------------------------------------------------------------------
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/loop/TStnAna.hh"

#include "Stntuple/obj/TStnGoodRunList.hh"

#include "Stntuple/scripts/init_geometry.C"
#include "Stntuple/scripts/setup_trigger_path.C" //    (const char* TriggerPath);
#include "Stntuple/scripts/parse_job_parameters.C" //  (TString& Parameters, StnAnaGlobals_t& Glob);

StnAnaGlobals_t         g;

void stnana (TString     Book   , 
	     TString     Dataset,
	     TString     Fileset       = "",
	     TString     File          = "",
	     TString     JobName       = "lumi()",      //
	     Int_t       NEvents       =  0)	 {    // 0: process all events
//-----------------------------------------------------------------------------
// step 1: make sure that all the needed scripts are loaded
// first read infrastructure scripts (global_vars.cc ... create_modules.C)
// it is assumed that 
// - there is a list of packages: pkg1, pkg2,pkg3 ...
// - each of the packages has a ana/scripts subdirectory;
// - a package 'pkg' has an ana/scripts/load_stnana_scripts_pkg.C script
//   which loads the jobs...
//-----------------------------------------------------------------------------
  char        macro[200], load_script[200], line[200], cmd[200];
  const char* pkg;

  TString test_release_dir = gEnv->GetValue("Stnana.TestReleaseDir",gSystem->Getenv("PWD"));

  printf("[stnana.C]: test_release_dir: %s\n",test_release_dir.Data());

  const char* script[] = { // loaded explicitly
    "init_geometry.C"       , "PWD",
    "parse_job_parameters.C", "PWD",
    "setup_trigger_path.C"  , "PWD",
    "init_photos.C"         , "STNTUPLE_MC_GEN",
    0 
  };

  TInterpreter::EErrorCode rc;

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i+=2) {
    const char* dir = gSystem->Getenv(script[i+1]);
    if (dir) {
      sprintf(macro,"%s/Stntuple/scripts/%s",dir,script[i]);
      if (! cint->IsLoaded(macro)) {
	printf("[stnana.C] locating %s\n",macro);
	cint->LoadMacro(macro,&rc);
	if (rc != 0) printf("stnana ERROR : failed to load %s\n",macro);
      }
    }
    else {
      printf("[stnana.C] WARNING: environment variable %s is not defined\n",script[i+1]);
    }
  }

  // if (! cint->IsLoaded(macro)) {
  //   cint->LoadMacro(macro);
  // }

  TString stnana_packages = gEnv->GetValue("Stnana.Package","");
  
  printf("[stnana.C] stnana packages: %s\n",stnana_packages.Data());

  TObjArray* list_of_packages = stnana_packages.Tokenize(" ");
  
  for (int i=0; i<list_of_packages->GetEntries(); i++) {
    pkg = ((TObjString*) list_of_packages->At(i))->String().Data();

    sprintf(macro,"%s/%s/ana/scripts/load_stnana_scripts_%s.C",test_release_dir.Data(),pkg,pkg);

    //    printf("[stnana.C] loading: %s\n",macro);

    if (! gInterpreter->IsLoaded(macro)) {
      gInterpreter->LoadMacro(macro);
      sprintf(load_script,"load_stnana_scripts_%s();",pkg);
      gInterpreter->ProcessLine(load_script);
    }
  }
//-----------------------------------------------------------------------------
// make sure necessary environment variables are defined consistently
//-----------------------------------------------------------------------------
  TString  dsid = Dataset;
  dsid.ReplaceAll('/','_');
  gSystem->Setenv("DSID"   ,dsid.Data());
  gSystem->Setenv("FILESET",Fileset);
//-----------------------------------------------------------------------------
// parse job options: /mc[=] /grl= /little /newcuts /output[=] /save[=] /debug= /pass=
//-----------------------------------------------------------------------------
  printf("[stnana.C] parsing command line\n");
  parse_job_parameters(JobName,g);

  if (g.dataset != 0) {
//-----------------------------------------------------------------------------
// deletion of the analysis loop has to be followed by deletion of all the modules
//-----------------------------------------------------------------------------
    delete g.x; 
    delete g.catalog;

    g.x            = 0;
    g.dataset      = 0;
    g.catalog      = 0;
    g.gGoodRunList = 0;
  }
  if (! g.dataset) {
//-----------------------------------------------------------------------------
//  initialize catalog and dataset, do it only once
//-----------------------------------------------------------------------------
    g.catalog = new TStnCatalog();
    g.dataset = new TStnDataset();

    if (Book == "pythia") {
//-----------------------------------------------------------------------------
// generator-level MC study - use Pythia, this approach can be extended  
// to work with other generators
//-----------------------------------------------------------------------------
      printf(">>> STNANA: PYTHIA initialization temporarily turned OFF\n");
      g.x = new TStnAna();
      // py  = TG3Pythia6::Instance();
      // m_gen = new TStnGeneratorModule();
      // m_gen->AddGenerator(py);
      // g.x->SetInputModule(m_gen);
    }
    else if (Book == "photos") {
//-----------------------------------------------------------------------------
// generator-level MC study - use PHOTOS
//-----------------------------------------------------------------------------
      printf("[stnana.C] PHOTOS initialization\n");
      g.x = new TStnAna();
      gInterpreter->ProcessLine("init_photos()");
    }
    else if (Book == "script") {
      g.x = new TStnAna();
    }
    else {
//-----------------------------------------------------------------------------
// analysis job runs on an exiting dataset. In case a good run list defines 
// a run range, only files from this run range get included
//-----------------------------------------------------------------------------
      if (Book == "") Book = "file";

/* DEBUG */     printf("[stnana.C]: Book=%s  Dataset=%s",Book.Data(),Dataset.Data());
/* DEBUG */     printf(" Fileset=%s  File=%s\n",Fileset.Data(),File.Data());
      /* DEBUG */  //   printf("[stnana.C]: g.MinRun=%i  g.MaxRun=%i\n",g.MinRun,g.MaxRun);    

      g.catalog->InitDataset(g.dataset,Book,Dataset,Fileset,File,g.MinRun,g.MaxRun);

      if (g.dataset->GetNFiles() <= 0) {
	printf(" empty dataset %s! exiting...\n",g.dataset->GetName());
	return;
      }

      /* DEBUG */     // printf("[stnana.C] done with the dataset initialization\n");
//-----------------------------------------------------------------------------
//  no matter what command line prevails
//-----------------------------------------------------------------------------
      if      (g.DoMc                 != -1) g.dataset->SetMcFlag(g.DoMc);
      else if (g.dataset->GetMcFlag() != -1) g.DoMc = g.dataset->GetMcFlag();
      else {
	g.DoMc = 0;
	g.dataset->SetMcFlag(0);
      }
      
      printf("[stnana.C] dataset MC_FLAG = %i\n",g.dataset->GetMcFlag());

      g.x = new TStnAna(g.dataset);
      g.x->SetPrintLevel(0);
    }
    
    // /* * DEBUG */    g.x->SetPrintLevel(100);
    g.x->GetInputModule()->SetPrintLevel(1);
//-----------------------------------------------------------------------------
// initialize good run list to be used
// so far assume that it is stored locally or on the disk pool
//-----------------------------------------------------------------------------
    if (g.GoodRunList != "NONE") {
      g.gGoodRunList = new TStnGoodRunList(g.GoodRunList.Data());
      g.x->SetGoodRunList(g.gGoodRunList);
    }
//-----------------------------------------------------------------------------
// initialize geometry, no action by default - Stntuple/scripts/init_geometry.C
// provides an empty routine
//-----------------------------------------------------------------------------
    const char* init_geometry = gEnv->GetValue("Stnana.InitGeometry",
					       "stntuple_init_geometry");
    sprintf(line,"%s();",init_geometry);
    printf("[stnana.C] executing:  %s\n",line);
    gInterpreter->ProcessLine(line);
    printf("[stnana.C] done with executing:  %s\n",line);
  }
  //  /* DEBUG */      return;
//-----------------------------------------------------------------------------
// setup trigger path check if required. When running on MC dataset need 
// to specify '/mc' in either batch or interactive job to make trigger 
// emulation to work
//-----------------------------------------------------------------------------
  printf("[stnana.C] starting trigger setup\n");    

  setup_trigger_path(g.L3TrigPath.Data());

  // /* DEBUG */     printf("[stnana.C] after setup_trigger_path\n");    
//-----------------------------------------------------------------------------
//  analyse definition of the requested job, handle debug mode
//-----------------------------------------------------------------------------
  // printf(" --- job_name = .%s. task = %s\n",g.JobName.Data(), g.Task.Data());
  
  int ind       = g.Task.Index("(");

  TString task  = g.Task(0,ind);
  TObjString* s = (TObjString*) g.ListOfTasks->FindObject(task.Data());
  if (s) {
//-----------------------------------------------------------------------------
// configure the job
//-----------------------------------------------------------------------------
    sprintf(cmd,"%s;",g.Task.Data());
    printf("cmd=%s\n",cmd);
    gInterpreter->ProcessLine(cmd,&rc);
    if (rc < 0) {
      printf("[stnana.C] called script returned -1, %s, bailing out ****** \n",task.Data());
      return;
    }
  }
  else {
    printf("[stnana.C] ERROR: unknown job : %s, bailing out ****** \n",task.Data());
    g.ListOfTasks->Print();
    return;
  }
  
  //  if (g.Debug != 0) debug(m_dbg);
//-----------------------------------------------------------------------------
//  output module if requested
//-----------------------------------------------------------------------------
  TStnOutputModule* om;

  if (g.OutputFileName != "") {
    printf(" ... writing output file %s\n",g.OutputFileName.Data());
    om = new TStnOutputModule(g.OutputFileName.Data());
    om->SetMaxFileSize(1800);
    g.x->SetOutputModule(om);
  }
//-----------------------------------------------------------------------------
// run the job, time its execution
//-----------------------------------------------------------------------------
  TStopwatch t;
  t.Start();
  if      (g.JobType == 0) g.x->Run(NEvents,g.MinRun,g.MaxRun);
  else if (g.JobType == 1) g.x->ProcessEventList(g.EventList);
  else if (g.JobType == 2) g.x->ProcessEventList(g.RunEventList);
  t.Stop();
  t.Print();
//-----------------------------------------------------------------------------
//  end of job: save histogram file, from now on now in Mode=2
//-----------------------------------------------------------------------------
  if (g.HistFileName != "") {
    printf(" ... saving histograms into %s\n",g.HistFileName.Data());
    g.x->SaveHist(g.HistFileName.Data(),2);
  }
//-----------------------------------------------------------------------------
// cleanup
//-----------------------------------------------------------------------------
  list_of_packages->Delete();
}
