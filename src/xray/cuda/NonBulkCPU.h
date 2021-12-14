#ifndef AMBER_NON_BULK_CPU_H
#define AMBER_NON_BULK_CPU_H

#include "xray/NonBulk.h"

namespace xray {
  class NonBulkCPU : public NonBulk {
  public:
    using NonBulk::NonBulk;

    ~NonBulkCPU() override = default;

    void calc_f_non_bulk(
      int n_atoms,
      const double* frac_xyz
    ) override;

  };
}
#endif //AMBER_NON_BULK_CPU_H
