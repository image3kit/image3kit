#pragma once

#ifdef PYBINDINGSTREAM
#include "globals.h"
#include <string>

#include <pybind11/pybind11.h>
namespace py = pybind11;
namespace pyd = pybind11::detail;

class __attribute__ ((visibility("default"))) instream {
private:
    py::dict            dic_;
    pyd::dict_iterator  ind_;
    bool                oky_ = true;
public:
    instream() {}
    instream(py::object obj):
        dic_(obj), ind_(dic_.begin()) {}

    template<typename T>
    instream& operator >> (T& t) {
        if (!oky_ || ind_ == dic_.end()) { oky_ = false; return *this; }
        t = ind_->second().template cast<T>();
        std::cerr<<t<<std::endl;
        ++ind_;
        return *this;
    }
    template<typename T>
    instream& operator >> (var3<T>& t) {
        if (!oky_ || ind_ == dic_.end()) { oky_ = false; return *this; }
        py::tuple ttt = ind_->second();
        t = {ttt[0].cast<T>(), ttt[1].cast<T>(), ttt[2].cast<T>()};
        std::cerr<<t<<std::endl;
        ++ind_;
        return *this;
    }
    operator bool() const { return oky_; }
    void str(const std::string& s) { std::cerr<<s<<std::endl; dic_["?"] = s; }
    char peek() const { std::cerr<<dic_.begin()->first().cast<std::string>()<<std::endl;  return '?';
        // return dic_.begin();
    }
    bool good() const { return oky_; }
};
#else
#include <sstream>
typedef std::stringstream instream;
#endif
