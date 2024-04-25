#ifndef TRACKSCAPE_UTILS_HPP
#define TRACKSCAPE_UTILS_HPP

namespace DAGGER {

template <class fT> class BasePropStorer {

public:
  fT prop = 0;

  BasePropStorer() { ; }
  BasePropStorer(fT prop) { this->prop = prop; }

  static BasePropStorer<fT> create() { return BasePropStorer<fT>(); }
  static BasePropStorer<fT> create(fT prop) { return BasePropStorer<fT>(prop); }
  template <class T, class V>
  static void mix(fT w1, BasePropStorer<T> &prop1, fT w2,
                  BasePropStorer<V> &prop2) {
    if (w1 + w2 == 0)
      return;
    prop1.prop = (w1 * prop1.prop + w2 * prop2.prop) / (w1 + w2);
  }
};

template <class fT, class stored_t> class VerticalStorer {

public:
  int nnodes = 0;
  fT dz = 1.;
  std::vector<fT> topcellz;
  std::vector<std::vector<stored_t>> pile;

  VerticalStorer() { ; }
  VerticalStorer(fT dz, int nnodes) {
    this->nnodes = nnodes;
    this->dz = dz;
    this->topcellz = std::vector<fT>(this->nnodes, 0.);
    this->pile = std::vector<std::vector<stored_t>>(this->nnodes,
                                                    std::vector<stored_t>());
  }

  void pile_up(int i, fT zadd, stored_t &prop) {
    if (zadd == 0)
      return;

    if (this->pile[i].size() == 0)
      this->pile[i].emplace_back(stored_t::create());

    while (zadd > 0) {
      fT nextspace = this->dz - this->topcellz[i];
      fT tzused = std::min(zadd, nextspace);
      zadd -= tzused;
      stored_t::mix(this->topcellz[i], this->pile[i].back(), tzused, prop);

      this->topcellz[i] += tzused;

      if (zadd > 0) {
        this->topcellz[i] = 0;
        this->pile[i].emplace_back(stored_t::create());
      }
    }
  }

  stored_t remove(int i, fT zrem) {
    auto removed = stored_t::create();
    if (zrem == 0)
      return removed;

    fT cumrem = 0;

    while (zrem > 0) {
      fT torem = std::min(zrem, this->topcellz[i]);
      this->topcellz[i] -= torem;
      stored_t::mix(cumrem, removed, torem, this->pile[i].back());
      cumrem += torem;
      if (this->topcellz[i] == 0) {
        this->pile[i].pop_back();
        this->topcellz[i] = this->dz;
        if (this->pile[i].size() == 0) {
          this->topcellz[i] = 0;
          this->pile[i].emplace_back(stored_t::create());
          return removed;
        }
      }
      zrem -= torem;
    }
    return removed;
  }
};

} // namespace DAGGER

#endif
