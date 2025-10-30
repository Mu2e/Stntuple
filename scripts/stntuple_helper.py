#!/usr/bin/env python
#
#
import os, re, string, subprocess
#------------------------------------------------------------------------------
class stntuple_helper:
    """stntuple_helper: class to build stntuple code"""
#   This appears to behave like c++ static member and is initialized at class defintion time.
    sourceroot =  os.path.abspath('.')

    def __init__(self, env, debug=False):
        self._debug  = debug
        self._env    = env;
        self._list_of_object_files = [];

        if (debug) : print("[stntuple_helper::__init__] Dir:"+self._env.Dir('.').abspath)

        self.dd      = re.search('[^/]*/[^/]*/[^/]*$',self._env.Dir('.').abspath).group(0)
        self.dirname = self.dd.split('/')[0];   # THIS
        self.subdir  = self.dd.split('/')[2];
        self.libname = self.dirname+'_'+self.subdir
        self.suffix  = ".hh" ;
        self.tmp_dir = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/tmp/'+self.subdir;
        self._debug  = debug
        if (debug) : 
            print("[stntuple_helper::__init__] -------------- building directory: "+self.dirname+'/'+self.subdir)
#
#   Accesors
#
    def base(self):
        return self.sourceroot

    def tmp_dira(self):
        return self.tmp_dir

    def debug(self):
        return self._debug

    def do_codegen(self,cc,script):
        self._env.StntupleCodegen(cc,script);
        obj = cc.replace(".cc",".o");
        self._env.SharedObject(obj,cc)
# don't need this any more
#        if (obj.find("_module") < 0):
#            self._list_of_object_files.append(obj);

    def handle_dictionaries(self,suffix = ".hh",skip_list = []):
        self.suffix = suffix ;
#------------------------------------------------------------------------------
# generate dictionaries
#------------------------------------------------------------------------------
        list_of_linkdef_files = self._env.Glob(self.subdir+'/dict/*_linkdef.h', strings=True)
        if (self._debug):
            print ("[stntuple_helper::handle_dictionaries] ["+self.dirname+"/"+self.subdir+"] handle_dictionaries: list_of_linkdef_files = ",
                   list_of_linkdef_files,' skip_list:',skip_list)
            
        list_of_dict_files    = []

        for f in list_of_linkdef_files:
            linkdef       = f.split('/');
            linkdef_fn    = linkdef[len(linkdef)-1];

            if (self._debug): 
                print ("[stntuple_helper::handle_dictionaries] linkdef_fn = ",linkdef_fn, '\n skip_list:',skip_list)

            if (not linkdef_fn in skip_list):
                clname        = linkdef_fn.replace("_linkdef.h","");
                include       = self.subdir+'/'+clname+self.suffix;
            
                dict          = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/tmp/'+self.subdir+'/'+clname+'_dict.cxx';
                list_of_dict_files.append(dict);

                if (self._debug):
                    print ("[stntuple_helper::handle_dictionaries] linkdef = ",linkdef)
                    print ("[stntuple_helper::handle_dictionaries] include = ",include)

                self._env.StntupleRootCint(dict,[f,include])
#------------------------------------------------------------------------------
# compile dictionaries
#------------------------------------------------------------------------------
        # list   = [];
        if (self._debug): 
            print ("[stntuple_helper::handle_dictionaries] list_of_dict_files = ",list_of_dict_files)

        for dict in list_of_dict_files:
            # print("----- dict : ",dict);
            obj_cxx_file = dict.replace(".cxx",".o");

            if (dict.find("_module_dict") < 0) : 
                # list.append(obj_cxx_file);
                self._list_of_object_files.append(obj_cxx_file);
        
            include  = dict.replace(".cxx",".h");
            self._env.SharedObject(obj_cxx_file,dict)


    def compile_fortran(self,list_of_f_files, skip_list = []):
        if (self._debug):
            print  (self._dirname+"[build_libs]: list_of_f_files:"+self._subdir,list_of_f_files);

        for f in list_of_f_files:
            if (not f in skip_list):
                o = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/tmp/'+self.subdir+'/'+f.split('.')[0]+'.o'
                self._list_of_object_files.append(o);
                self._env.SharedObject(o,f)


    def build_libs(self,list_of_cc_files, skip_list = [],libs = []):
        if (self._debug):
            print ("["+self.dirname+".build_libs]: list_of_cc_files:"+self.subdir,list_of_cc_files)

        for cc in list_of_cc_files:
            if (not cc in skip_list):
                if (self._debug):
                    print ("stntuple_helper::build_libs: .cc file: "+cc)
                    
                o = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/tmp/'+self.subdir+'/'+cc.split('.')[0]+'.o'
                self._list_of_object_files.append(o);
                self._env.SharedObject(o,cc)
        #------------------------------------------------------------------------------
        # need to keep MU2E_SATELLITE_RELEASE for a while for building in a satellite release
        #------------------------------------------------------------------------------
        dir = os.environ.get('MUSE_BUILD_DIR')+'/'+self.dirname;
        lib_name = dir+'/lib/'+self.libname+'.so';

        if (self._debug):
            print ("list_of_object_files:",self._list_of_object_files)

        self._env.SharedLibrary(lib_name,self._list_of_object_files,LIBS = [libs])


    def build_modules(self,list_of_module_files, skip_list, skip_dict_list = [], libs = []):
        if (self._debug):
            print ("[Stntuple.build_modules]: list_of_module_files in Stntuple/"+self.subdir," : ",list_of_module_files)

        for module in list_of_module_files:
            if (not module in skip_list):
                if (self._debug):
                    print ("module file: "+module)
                    
                o = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/tmp/'+self.subdir+'/'+module.split('.')[0]+'.o'
                self._env.SharedObject(o,module)
                olist = [o];

                # check for the presence of a dictionary
                dict_bn = module.split('.')[0]+'_linkdef.h';
                if (not dict_bn in  skip_dict_list):
                    dict = self.dirname+'/'+self.subdir+'/'+self.subdir+'/dict/'+dict_bn;

                    if (os.path.isfile(os.getenv("PWD")+'/'+dict)):
                        odict = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/tmp/'+self.subdir+'/'+module.split('.')[0]+'_dict.o'
                        olist.append(odict)
                    # else: 
                        #    print(" no ",dict ," found in",os.getenv('PWD'));

                mname = os.path.basename(module).split('.')[0];
                lib   = os.getenv("MUSE_BUILD_DIR")+'/'+self.dirname+'/lib/libmu2e_'+self.dirname+'_'+mname+'.so';

                if (self._debug):
                    print ("o: "+o, "lib:"+lib)
                    print ("pwd: "+os.getenv("PWD"),"module dict: "+dict);
                    print ("olist : ",olist);

                self._env.SharedLibrary(lib,olist,LIBS =  ['libStntuple_mod.so',libs]);

    def print():
        print("stntuple_helper::print :  emoe");
#------------------------------------------------------------------------------ END
