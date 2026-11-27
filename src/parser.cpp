/*------------------------------------------------------------------------*/
/*! \file parser.cpp
    \brief contains functions necessary to parse the AIG

  Part of AMulet2 : AIG Multiplier Verification Tool.
  Copyright(C) 2020, 2021 Daniela Kaufmann, Johannes Kepler University Linz
*/
/*------------------------------------------------------------------------*/
#include "parser.h"

#include <stdint.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>
/*------------------------------------------------------------------------*/
// ERROR CODES:
static int err_parsing       = 20; // general parsing error
static int err_latches       = 21; // cannot handle latches
static int err_no_inputs     = 22; // no inputs
static int err_odd_inputs    = 23; // odd inputs
static int err_no_outputs    = 24; // no outputs
static int err_wrong_outputs = 25; // wrong number of outputs
static int err_rect_mult     = 26; // operands of different width
/*------------------------------------------------------------------------*/
// Operand shape check by two-hot simulation.
//
// Only the number of inputs NN is known, and the operands are assumed to have
// NN/2 bits each. To detect multipliers with operands of different width
// (e.g. 16x12, which also has 28 inputs and 28 outputs) we simulate the AIG
// with exactly two inputs x, y set to 1. For a multiplier a * b the result is
//   0                       if x and y belong to the same operand
//   +-2^(w(x)+w(y))         if they belong to different operands,
// where w is the bit weight (negative only for signed multipliers, if exactly
// one of x, y is a sign bit). This identifies both operands and all weights,
// independently of the multiplier architecture (array, tree, booth, ...).
// Both signs are accepted in either mode, since only the shape is checked.
/*------------------------------------------------------------------------*/
namespace {

typedef std::pair<unsigned, unsigned> InputPair;

const int TWO_HOT_ZERO      = -1;  // product is zero
const int TWO_HOT_MALFORMED = -2;  // not of the form +-2^k

class TwoHotSimulator {
  std::vector<unsigned> order;  // indices of AND gates in topological order
  std::vector<uint64_t> val;    // simulation value of each variable

  uint64_t lit_val(unsigned lit) const {
    uint64_t v = val[aiger_lit2var(lit)];
    return aiger_sign(lit) ? ~v : v;
  }

 public:
  TwoHotSimulator() : val(M, 0) {
    // AND gates in ASCII AIGER files do not need to be sorted
    std::vector<int> and_of(M, -1);
    for (unsigned i = 0; i < get_model_num_ands(); i++)
      and_of[aiger_lit2var(get_model_and(i)->lhs)] = i;

    std::vector<char> mark(M, 0);  // 0 = new, 1 = open, 2 = done
    std::vector<unsigned> stack;
    for (unsigned o = 0; o < NN; o++) {
      stack.push_back(aiger_lit2var(slit(o)));
      while (!stack.empty()) {
        unsigned v = stack.back();
        int idx = and_of[v];
        if (mark[v] == 2 || idx < 0) {
          mark[v] = 2;
          stack.pop_back();
        } else if (mark[v] == 1) {
          mark[v] = 2;
          order.push_back(idx);
          stack.pop_back();
        } else {
          mark[v] = 1;
          const aiger_and * a = get_model_and(idx);
          unsigned c0 = aiger_lit2var(a->rhs0), c1 = aiger_lit2var(a->rhs1);
          if (!mark[c0]) stack.push_back(c0);
          if (!mark[c1]) stack.push_back(c1);
        }
      }
    }
  }

  // For each pair, sets both inputs to 1 (all others to 0) and classifies
  // the output vector: TWO_HOT_ZERO, TWO_HOT_MALFORMED, or k for +-2^k.
  std::vector<int> run(const std::vector<InputPair> & pairs) {
    std::vector<int> res;
    res.reserve(pairs.size());
    std::vector<uint64_t> out(NN);

    for (size_t base = 0; base < pairs.size(); base += 64) {
      size_t lanes = std::min<size_t>(64, pairs.size() - base);

      for (unsigned i = 0; i < NN; i++)
        val[aiger_lit2var(get_model_inputs_lit(i))] = 0;
      for (size_t l = 0; l < lanes; l++) {
        const InputPair & p = pairs[base + l];
        val[aiger_lit2var(get_model_inputs_lit(p.first))]  |= 1ull << l;
        val[aiger_lit2var(get_model_inputs_lit(p.second))] |= 1ull << l;
      }
      for (unsigned idx : order) {
        const aiger_and * a = get_model_and(idx);
        val[aiger_lit2var(a->lhs)] = lit_val(a->rhs0) & lit_val(a->rhs1);
      }
      for (unsigned o = 0; o < NN; o++) out[o] = lit_val(slit(o));

      for (size_t l = 0; l < lanes; l++) {
        unsigned low = 0;
        while (low < NN && !((out[low] >> l) & 1)) low++;
        if (low == NN) { res.push_back(TWO_HOT_ZERO); continue; }

        // above 'low' all bits are 0 (+2^low) or all 1 (-2^low)
        bool all0 = true, all1 = true;
        for (unsigned o = low + 1; o < NN; o++) {
          if ((out[o] >> l) & 1) all0 = false;
          else all1 = false;
        }
        if (all0 || all1) res.push_back(low);
        else res.push_back(TWO_HOT_MALFORMED);
      }
    }
    return res;
  }
};

/*------------------------------------------------------------------------*/
// Compact list of input indices, e.g. "0..15" or "0,2,...,30"
std::string index_list(const std::vector<unsigned> & idx) {
  std::string s;
  if (idx.empty()) return s;
  if (idx.size() > 2) {
    unsigned d = idx[1] - idx[0];
    bool progression = true;
    for (size_t i = 2; i < idx.size(); i++)
      if (idx[i] - idx[i-1] != d) progression = false;
    if (progression) {
      if (d == 1)
        return std::to_string(idx[0]) + ".." + std::to_string(idx.back());
      return std::to_string(idx[0]) + "," + std::to_string(idx[1]) +
             ",...," + std::to_string(idx.back());
    }
  }
  for (size_t i = 0; i < idx.size(); i++) {
    if (i) s += ",";
    s += std::to_string(idx[i]);
  }
  return s;
}

/*------------------------------------------------------------------------*/
// Determines both operands by two-hot simulation. Returns false if the
// simulation results do not match any multiplier (e.g. a faulty circuit),
// in which case the default assumption of NN/2 x NN/2 is kept.
// For operands of different width, all pairs (or a fixed sample) are
// simulated in addition, and the first pair whose result does not match
// is stored in 'mismatch' (otherwise it is set to (0, 0)).
bool infer_operands(std::vector<unsigned> * op, InputPair * mismatch) {
  *mismatch = InputPair(0, 0);
  TwoHotSimulator sim;
  std::vector<int> side(NN, 0), w(NN, -1);
  std::vector<InputPair> pairs;

  // 1) input 0 together with every other input:
  //    nonzero result <=> the other input belongs to the other operand
  for (unsigned y = 1; y < NN; y++) pairs.push_back(InputPair(0, y));
  std::vector<int> res = sim.run(pairs);

  unsigned y0 = 0;
  int w0 = -1;
  for (unsigned y = 1; y < NN; y++) {
    int r = res[y-1];
    if (r == TWO_HOT_MALFORMED) return false;
    if (r == TWO_HOT_ZERO) continue;
    side[y] = 1;
    if (w0 < 0 || r < w0) w0 = r, y0 = y;
  }
  if (!y0) return false;
  w[0] = w0;  // the lowest result is w(input 0) + 0
  for (unsigned y = 1; y < NN; y++)
    if (side[y]) w[y] = res[y-1] - w0;

  // 2) y0 (weight 0) together with every other input of input 0's operand
  pairs.clear();
  for (unsigned x = 1; x < NN; x++)
    if (!side[x]) pairs.push_back(InputPair(x, y0));
  res = sim.run(pairs);
  for (size_t i = 0; i < pairs.size(); i++) {
    if (res[i] < 0) return false;
    w[pairs[i].first] = res[i];
  }

  // 3) each operand has exactly the weights 0, ..., width-1
  op[0].clear(), op[1].clear();
  for (unsigned i = 0; i < NN; i++) op[side[i]].push_back(i);
  for (int s = 0; s < 2; s++) {
    std::vector<char> seen(op[s].size(), 0);
    for (unsigned i : op[s]) {
      if (w[i] < 0 || static_cast<size_t>(w[i]) >= op[s].size()) return false;
      if (seen[w[i]]) return false;
      seen[w[i]] = 1;
    }
  }
  if (op[0].size() == op[1].size()) return true;

  // 4) different widths: check all pairs (or a fixed sample of them),
  //    to tell a correct rectangular multiplier from a faulty one
  const size_t max_pairs = 64 * 256;
  size_t num_pairs = static_cast<size_t>(NN) * (NN-1) / 2;
  pairs.clear();
  if (num_pairs <= max_pairs) {
    for (unsigned x = 0; x < NN; x++)
      for (unsigned y = x + 1; y < NN; y++) pairs.push_back(InputPair(x, y));
  } else {
    uint64_t seed = 88172645463325252ull;  // xorshift64
    while (pairs.size() < max_pairs) {
      seed ^= seed << 13, seed ^= seed >> 7, seed ^= seed << 17;
      unsigned x = seed % NN, y = (seed >> 32) % NN;
      if (x != y) pairs.push_back(InputPair(x, y));
    }
  }
  res = sim.run(pairs);
  for (size_t i = 0; i < pairs.size(); i++) {
    unsigned x = pairs[i].first, y = pairs[i].second;
    int expected = side[x] == side[y] ? TWO_HOT_ZERO : w[x] + w[y];
    if (res[i] != expected) {
      *mismatch = InputPair(std::min(x, y), std::max(x, y));
      break;
    }
  }
  return true;
}

}  // namespace

/*------------------------------------------------------------------------*/
// Rejects multipliers whose operands have different widths.
static void check_operand_widths() {
  std::vector<unsigned> op[2];
  InputPair mismatch;
  if (!infer_operands(op, &mismatch)) {
    if (verbose >= 2)
      msg("operand shape check inconclusive, assuming %ux%u multiplier",
        NN/2, NN/2);
    return;
  }

  if (op[0].size() < op[1].size()) std::swap(op[0], op[1]);
  unsigned na = op[0].size(), nb = op[1].size();
  if (na == nb) {
    if (verbose >= 2) msg("operand shape check: %ux%u multiplier", na, nb);
    return;
  }

  std::string faulty;
  if (mismatch.first != mismatch.second)
    faulty = "\n*** [amulet2] the circuit is also not a correct multiplier "
             "of this shape (simulation with only inputs " +
             std::to_string(mismatch.first) + " and " +
             std::to_string(mismatch.second) + " set to 1 does not match)";

  die(err_rect_mult,
    "detected a %ux%u multiplier (operands: inputs %s and inputs %s)%s\n"
    "*** [amulet2] only multipliers with operands of equal width are supported",
    na, nb, index_list(op[0]).c_str(), index_list(op[1]).c_str(),
    faulty.c_str());
}

/*------------------------------------------------------------------------*/

bool match_and(unsigned lhs, unsigned rhs0, unsigned rhs1) {
  if (lhs == aiger_false) return 0;
  if (aiger_sign(lhs)) return 0;
  assert(lhs != aiger_true);

  aiger_and * and1 = is_model_and(lhs);
  if (!and1) return 0;
  if (and1->rhs0 == rhs0 && and1->rhs1 == rhs1) return 1;
  if (and1->rhs0 == rhs1 && and1->rhs1 == rhs0) return 1;
  return 0;
}

/*------------------------------------------------------------------------*/
void determine_input_order() {
  unsigned s0 = 0, sl = NN-1;
  if (match_and(
        slit(0),
        get_model_inputs_lit(0),
        get_model_inputs_lit(1))) {
    a0 = 0, al = NN-2, ainc = 2;
    b0 = 1, bl = NN-1, binc = 2;
    msg("assuming ordering as in BTOR generated benchmarks");
  } else {
    a0 = 0, al = NN/2-1,   ainc = 1;
    b0 = NN/2, bl = NN-1, binc = 1;
    msg("assuming ordering as in the ABC generated or AOKI benchmarks");
  }
  if (verbose >= 2) {
    if (NN == 2) {
      msg("a[0] = input[%d]", a0);
      msg("b[0] = input[%d]", b0);
      msg("s[0] = output[%d]", s0);
    } else if (NN == 4) {
      msg("(a[0], a[1]) =(input[%d], input[%d])", a0, al);
      msg("(b[0], b[1]) =(input[%d], input[%d])", b0, bl);
      msg("(s[0], ..., s[3]) =(output[%d], ..., output[%d])", s0, sl);
    } else if (NN == 6) {
      msg("(a[0], a[1], a[2]) =(input[%d], input[%d], input[%d])",
        a0, a0 + ainc, al);
      msg("(b[0], b[1], b[2]) =(input[%d], input[%d], input[%d])",
        b0, b0 + binc, bl);
      msg("(s[0], ..., s[5]) =(output[%d], ..., output[%d])", s0, sl);
    } else {
      msg("(a[0], a[1], ..., a[%d]) =(input[%d], input[%d], ..., input[%d])",
        NN/2-1, a0, a0 + ainc, al);
      msg("(b[0], b[1], ..., b[%d]) =(input[%d], input[%d], ..., input[%d])",
        NN/2-1, b0, b0 + binc, bl);
      msg("(s[0], ..., s[%d]) =(output[%d], ..., output[%d])",
        NN-1, s0, sl);
    }
  }
}

/*------------------------------------------------------------------------*/

void init_aiger_with_checks() {
  if (get_model_num_latches()) die(err_latches, "can not handle latches");
  if (!get_model_num_inputs()) die(err_no_inputs, "no inputs");
  if ((get_model_num_inputs() & 1)) die(err_odd_inputs, "odd number of inputs");
  if (!get_model_num_outputs()) die(err_no_outputs, "no outputs");
  if (get_model_num_outputs() == get_model_num_inputs()) {
    M = get_model_maxvar() + 1;
    NN = get_model_num_outputs();
  }
  else  die(err_wrong_outputs, "expected %u but got %u outputs",
      get_model_num_inputs(), get_model_num_outputs());


  msg("MILOA %u %u %u %u %u",
    get_model_maxvar(),
    get_model_num_inputs(),
    get_model_num_latches(),
    get_model_num_outputs(),
    get_model_num_ands());

  check_operand_widths();
  determine_input_order();
}

/*------------------------------------------------------------------------*/

void parse_aig(const char * input_name) {
  init_aig_parsing();

  msg("reading '%s'", input_name);
  const char * err = aiger_open_and_read_to_model(input_name);
  if (err) die(err_parsing, "error parsing '%s': %s", input_name, err);

  init_aiger_with_checks();
}
