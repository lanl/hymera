#pragma once

#include "Tasks.h"

TaskStatus ScaleView(View3 view, const Real scl);
TaskStatus ReduceView(View3 view);
TaskStatus AccumulateView(DualView3 acc_view, View3 view);
TaskStatus ZeroView(View3 view);
