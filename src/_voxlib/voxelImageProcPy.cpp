#include <pybind11/pybind11.h>

#include <string_view>
#include <iostream>
#include <string>

// pk_default_print
namespace py = pybind11;
inline py::arg operator""_a (const char* c, size_t) { return py::arg(c); }
inline py::str operator""_s (const char* c, size_t) { return py::str(c); }

using namespace std;
int main(int argc, char** argv) {

  try {
        assert(argc > 1 && "Usage:\n pocketpy -c 'commands' \nor\n pocketpy file.py");

    py::initialize(); // py_initialize();
    //py_sys_setargv(argc, argv);


    //auto sim = py::module::import("smod");
    //auto result = sim.attr("srun")("echo 123sdsXXXsddsd").cast<int>();
    //assert(result == 0);
        // std::cout << smod.attr("__version__").cast<std::string>() << std::endl;
    auto locals = py::dict();
    auto glob = py::globals();
        // py::setattr(locals,"smod", smod);
        // py::setattr(locals, "smod", py_getdict(smod.ptr(), py_name("smod")));

        // locals["smod"_s] = smod;


    //py::print(py::globals());
    //int ret = py::cast<int>(sir("echo test call SiR from python wrapped in SiM!"));
    if (argc > 2) {
      try {
        // py::exec("from smod import voxelImageI", glob);
        py::exec(argv[2], glob, locals);//py_exec(ins.str().c_str(), glob);
        //py::module_ sys = py::module_::import("sys");
        //py::print(sys.attr("path"));
      }
      catch (std::exception& e)  {
        cout << "Error in py.exec: " << e.what()<<endl;
      }
      catch (...)  {
        cout << "Error in py.exec "<<endl;
      }
    }
        else {
            std::cout << "Error, usage:\n pocketpy -c 'commands' \n";
            return 1;
        }
  }
  catch (std::exception& e) {
    cout << "Error in py: " << e.what()<<endl;
  }
  catch (...) { cout << "Error in py: "<<endl; }

  int code = py_checkexc(false) ? 1 : 0;
  py::finalize(true); //py_finalize();
  return code;
}
