#pragma once

// Global #defines

#ifdef _MSC_VER
#   define CAP_FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__GNUG__) || defined(__clang__)
#   define CAP_FORCEINLINE __attribute((always_inline))
#else 
#   define CAP_FORCEINLINE inline
#endif

#ifdef NDEBUG
#   define MAX_SPEED
#endif
