#! /usr/bin/env python

from waflib.Configure import conf

@conf
def check_with(conf, check, what, *args, **kwargs):
    """
    Perform `check`, also looking at --with-X commandline option and
    and X_HOME environment variable
    """
    import os
    
    with_dir = getattr(conf.options, "with_" + what, None)
    env_dir = os.environ.get(what.upper() + "_HOME", None)
    paths = [with_dir, env_dir] + kwargs.pop("extra_paths", [])
    
    what = what.upper()
    kwargs["uselib_store"] = kwargs.get("uselib_store", what)
    kwargs["use"] = kwargs.get("use", []) + [kwargs["uselib_store"]]
    
    for path in [p for p in paths if p]:
        if conf.find_at(check, what, path, **kwargs):
            return
            
    check(**kwargs)
    
@conf
def find_at(conf, check, what, where, **kwargs):
    from os.path import exists, join as pjoin
    if not exists(where):
        return False
    try:
        conf.env.stash()
        conf.env[what + "_HOME"] = where
        conf.env.append_value('PATH',  pjoin(where, "bin"))
        conf.env.append_value('RPATH', pjoin(where, "lib"))
        conf.env.append_value('PKG_CONFIG_PATH', pjoin(where, "lib/pkgconfig"))
        conf.parse_flags("-I{0}/include -L{0}/lib".format(where), 
                         uselib=kwargs["uselib_store"])
        check(**kwargs)
        return True
    except conf.errors.ConfigurationError:
        conf.end_msg("failed", color="YELLOW")
        conf.env.revert()
        return False
        
def options(opt):
    opt.load('compiler_cxx')
    opt.add_option('--with-bat', default=None,
        help="Also look for BAT at the given path")
    
def configure(conf):
    conf.load('compiler_cxx')
    
    conf.check_with(conf.check_cfg, "cern_root_system",
                    path="root-config", args="--cflags --libs", package="")
                    
    conf.check_cxx(lib="Minuit", use=["CERN_ROOT_SYSTEM"],
                   uselib_store="CERN_ROOT_SYSTEM")
    
    conf.check_with(conf.check_cxx, "bat", 
                    stlib="BAT",
                    header="BAT/BCModel.h", 
                    use=["CERN_ROOT_SYSTEM"])
    
    zpdir = conf.path.find_node("ZPrimeBATTutorial/limits")
    conf.parse_flags("-I{0}".format(zpdir.abspath()), uselib="MTFA_INC")
    conf.check_with(conf.check_cxx, "mtfa",
                    header_name="BCMTFAnalysisFacility.h",
                    use=["MTFA_INC"])
    
    conf.to_log("Final environment:")
    conf.to_log(conf.env)
    
def build(bld):
    
    bld.stlib(source=bld.path.ant_glob("ZPrimeBATTutorial/limits/BC*.cxx"),
              target="MTFA", use=["BAT", "CERN_ROOT_SYSTEM", "MTFA_INC"])

    bld.program(
        source=bld.path.ant_glob("*.cxx"),
        target="main",
        use=["BAT", "CERN_ROOT_SYSTEM", "MTFA"])
