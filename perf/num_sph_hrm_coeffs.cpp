#include <cstddef>
#include <iostream>
#include <chrono>
#include <typeinfo>

template<typename size_type>
size_type num_sph_hrm_coeffs(size_type const l_max,size_type const m_max) {
    return (m_max+size_type(1))*(l_max+size_type(2))/size_type(2)
        +(l_max+size_type(1))*(l_max-m_max);
}

template<typename size_type>
void estimate_elapsed_time() {
    size_type const l_max(4096);
    size_type const m_max(l_max);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    size_t const repeat_times(100000000);
    for(size_t i=0;i<repeat_times;++i){
        num_sph_hrm_coeffs(l_max,m_max);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"Elapsed time for "<< typeid(l_max).name()
        << " " <<elapsed_seconds.count()/double(repeat_times)<<std::endl;
}

int main(){
    estimate_elapsed_time<size_t>();
    estimate_elapsed_time<ptrdiff_t>();
    estimate_elapsed_time<unsigned long>();
    estimate_elapsed_time<unsigned long long>();
    return 0;
}