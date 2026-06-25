#pragma once
#include "fid.h"
#include "FieldComponents.h"

KOKKOS_INLINE_FUNCTION
constexpr int component_offset(field_id id, int time) {
  switch (id) {
    case field_id::B:  return (time == 0) ? FieldComponents::B : FieldComponents::Bt;
    case field_id::J:  return (time == 0) ? FieldComponents::J : FieldComponents::Jt;
    case field_id::E:  return (time == 0) ? FieldComponents::E : FieldComponents::Et;
    default:           return -1;
  }
}


