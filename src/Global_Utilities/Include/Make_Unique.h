#ifndef __MAKE_UNIQUE_H__
#define __MAKE_UNIQUE_H__

#include <memory>

namespace fdaPDE{

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique_time(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

}
#endif
