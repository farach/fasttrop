/*
 * Operating-system helpers for the fasttrop core, kept in their own
 * translation unit so <windows.h> never meets the R headers.
 *
 * - trop_wtime(): high-resolution monotonic wall clock. On Windows,
 *   omp_get_wtime() from MinGW's libgomp has only millisecond resolution, so
 *   QueryPerformanceCounter() is used instead.
 * - trop_qos_enter()/trop_qos_leave(): on Windows, opt the calling thread out
 *   of execution-speed power throttling (EcoQoS) while a native kernel runs,
 *   then return control of the thread to the system. Windows 11 throttles
 *   long-running background processes (such as an R session doing a long
 *   computation without a foreground window), which slowed identical native
 *   work by about 1.3x in a measured R session. Calls nest per thread. They do
 *   nothing on other platforms or when the user allows throttling.
 */
#ifdef _WIN32
#if !defined(_WIN32_WINNT) || _WIN32_WINNT < 0x0602
#undef _WIN32_WINNT
#define _WIN32_WINNT 0x0602
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <time.h>
#endif

#include "trop_core.h"

double trop_wtime(void) {
#ifdef _WIN32
  LARGE_INTEGER f, c;
  QueryPerformanceFrequency(&f);
  QueryPerformanceCounter(&c);
  return (double) c.QuadPart / (double) f.QuadPart;
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
#endif
}

#ifdef _WIN32
/* Not declared by every MinGW-w64 release; layout from the Windows SDK. */
typedef struct {
  ULONG Version;
  ULONG ControlMask;
  ULONG StateMask;
} trop_thread_throttling_state;
#define TROP_THROTTLING_VERSION 1UL
#define TROP_THROTTLING_EXECUTION_SPEED 0x1UL
#define TROP_THREAD_POWER_THROTTLING ((THREAD_INFORMATION_CLASS) 3) /* ThreadPowerThrottling */

static _Thread_local int qos_depth = 0;

static void set_throttling(ULONG control_mask) {
  trop_thread_throttling_state s;
  s.Version = TROP_THROTTLING_VERSION;
  s.ControlMask = control_mask; /* 0 hands the decision back to the system */
  s.StateMask = 0;              /* with ControlMask set: do not throttle */
  SetThreadInformation(GetCurrentThread(), TROP_THREAD_POWER_THROTTLING, &s, sizeof s);
}
#endif

void trop_qos_enter(const trop_opts *o) {
#ifdef _WIN32
  if (o->throttle) return;
  if (qos_depth++ == 0) set_throttling(TROP_THROTTLING_EXECUTION_SPEED);
#else
  (void) o;
#endif
}

void trop_qos_leave(const trop_opts *o) {
#ifdef _WIN32
  if (o->throttle) return;
  if (qos_depth > 0 && --qos_depth == 0) set_throttling(0UL);
#else
  (void) o;
#endif
}
