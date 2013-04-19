#! /usr/bin/env python

from waflib.Configure import conf
        
def options(opt):
    opt.load('compiler_cxx')
    opt.load('compiler_magic check_with', tooldir="common/waf")
    
    opt.add_option('--with-bat', default=None,
        help="Also look for BAT at the given path")
        
    opt.add_option('--with-cern-root-system', default=None,
        help="Look for ROOT at the given path")
        
    opt.add_option('--fully-static', action="store_true")
    
def configure(conf):
    conf.load('compiler_cxx')
    conf.load('compiler_magic check_with', tooldir="common/waf")
    
    if conf.options.fully_static:
        conf.check_with(conf.check_cxx, "cern_root_system", 
                        stlib=["Root", "lzma", "pcre", "z", "m", "freetype"]) #, "crypt", "ssl", "krb5", "k5crypto", "krb5support", "crypto", "com_err", "resolv", "selinux", "sepol"])
                        
        conf.env.append_value("CXXFLAGS_CERN_ROOT_SYSTEM", ["-pthread"])
        conf.env.append_value("LINKFLAGS_CERN_ROOT_SYSTEM", ["-pthread"])
        conf.env.append_value("LIB_CERN_ROOT_SYSTEM", ["dl"])
        #conf.env.append_value("LIB_CERN_ROOT_SYSTEM", ["dl", "freetype"])
        
    else:
        conf.check_with(conf.check_cfg, "cern_root_system",
                        path="root-config", args="--cflags --libs", package="")
        conf.check_cxx(lib="Minuit", use=["CERN_ROOT_SYSTEM"],
                       uselib_store="CERN_ROOT_SYSTEM")
        conf.env.RPATH_CERN_ROOT_SYSTEM = conf.env.LIBPATH_CERN_ROOT_SYSTEM
        
    conf.env.append_value("LINKFLAGS", ["-Wl,--as-needed"])
    
    conf.check_with(conf.check_cxx, "bat", 
                    stlib="BAT",
                    header="BAT/BCModel.h", 
                    use=["CERN_ROOT_SYSTEM"])
    
    zpdir = conf.path.find_node("ZPrimeBATTutorial/limits")
    conf.parse_flags("-I%s" % zpdir.abspath(), uselib="MTFA_INC")
    conf.check_with(conf.check_cxx, "mtfa",
                    header_name="BCMTFAnalysisFacility.h",
                    use=["MTFA_INC"])
    
    conf.env.append_value("CXXFLAGS", ["-g", "-O2", "-Werror", "-Wall", "-pthread"])
    conf.env.append_value("LINKFLAGS", ["-pthread"])
    
    conf.env.append_value("STLIB", ["m", "c"])
    
    if conf.options.fully_static:
        #conf.env.SHLIB_MARKER = '-Wl,-Bstatic'
        conf.env.append_value("LINKFLAGS", ["-static-libgcc", "-static-libstdc++"])
        pass
    
    conf.to_log("Final environment:")
    conf.to_log(conf.env)
    
def build(bld):
    
    bld(features="cxx", source=bld.path.ant_glob("ZPrimeBATTutorial/limits/BC*.cxx"),
              target="MTFA", use=["BAT", "CERN_ROOT_SYSTEM", "MTFA_INC"])

    bld.program(
        source=bld.path.ant_glob("*.cxx"),
        target="diphoton2012-limitsetting",
        use=["BAT", "CERN_ROOT_SYSTEM", "MTFA"])
