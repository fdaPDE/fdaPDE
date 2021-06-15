#ifndef __MAKE_SHARED_H__
#define __MAKE_SHARED_H__

#include <memory>
namespace fdaPDE{
template<typename T, typename... Args>
std::shared_ptr<T> make_shared(Args&&... args)
{
    return std::shared_ptr<T>(new T(std::forward<Args>(args)...));
}
};

#endif
