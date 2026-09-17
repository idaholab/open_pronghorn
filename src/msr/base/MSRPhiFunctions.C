#include "MSRPhiFunctions.h"

#include "MooseError.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>

namespace MSR
{
namespace
{

DenseMatrix
identityMatrix(const std::size_t n)
{
  DenseMatrix I(n, std::vector<Real>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i)
    I[i][i] = 1.0;
  return I;
}

DenseMatrix
matMul(const DenseMatrix & A, const DenseMatrix & B)
{
  const std::size_t n = A.size();
  if (B.size() != n)
    mooseError("MSRPhiFunctions: matrix size mismatch.");

  DenseMatrix C(n, std::vector<Real>(n, 0.0));

  for (std::size_t i = 0; i < n; ++i)
  {
    if (A[i].size() != n || B[i].size() != n)
      mooseError("MSRPhiFunctions: matrix must be square.");

    for (std::size_t k = 0; k < n; ++k)
    {
      const Real aik = A[i][k];
      if (aik == 0.0)
        continue;

      for (std::size_t j = 0; j < n; ++j)
        C[i][j] += aik * B[k][j];
    }
  }

  return C;
}

Real
matrixInfinityNorm(const DenseMatrix & A)
{
  Real norm = 0.0;

  for (const auto & row : A)
  {
    Real row_sum = 0.0;
    for (const auto value : row)
      row_sum += std::abs(value);

    norm = std::max(norm, row_sum);
  }

  return norm;
}

DenseMatrix
matrixExponential(const DenseMatrix & A)
{
  const std::size_t n = A.size();

  for (const auto & row : A)
    if (row.size() != n)
      mooseError("MSRPhiFunctions: matrix must be square.");

  if (n == 0)
    return {};

  const Real norm = matrixInfinityNorm(A);

  unsigned int squarings = 0;
  if (norm > 0.5)
    squarings = static_cast<unsigned int>(std::ceil(std::log2(norm / 0.5)));

  const Real scale = std::ldexp(1.0, -static_cast<int>(squarings));

  DenseMatrix B = A;
  for (auto & row : B)
    for (auto & value : row)
      value *= scale;

  DenseMatrix result = identityMatrix(n);
  DenseMatrix term = identityMatrix(n);

  const Real tolerance = 1e-15;
  const unsigned int max_terms = 200;

  bool converged = false;

  for (unsigned int j = 1; j <= max_terms; ++j)
  {
    term = matMul(term, B);

    const Real inv_j = 1.0 / static_cast<Real>(j);
    for (auto & row : term)
      for (auto & value : row)
        value *= inv_j;

    Real term_norm = 0.0;
    Real result_norm = 0.0;

    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t k = 0; k < n; ++k)
      {
        result[i][k] += term[i][k];
        term_norm = std::max(term_norm, std::abs(term[i][k]));
        result_norm = std::max(result_norm, std::abs(result[i][k]));
      }

    if (term_norm <= tolerance * std::max(1.0, result_norm))
    {
      converged = true;
      break;
    }
  }

  if (!converged)
    mooseError("MSRPhiFunctions: scaled matrix exponential Taylor series did not converge.");

  for (unsigned int i = 0; i < squarings; ++i)
    result = matMul(result, result);

  return result;
}


struct SmallExpStats
{
  unsigned int terms = 0;
  unsigned int squarings = 0;
};

constexpr std::size_t MAX_SMALL_DIM = 16;
using SmallMatrix = std::array<Real, MAX_SMALL_DIM * MAX_SMALL_DIM>;
using SmallVector = std::array<Real, MAX_SMALL_DIM>;

inline std::size_t
smallIndex(const std::size_t i, const std::size_t j)
{
  return i * MAX_SMALL_DIM + j;
}

void
smallIdentity(SmallMatrix & A, const std::size_t n)
{
  A.fill(0.0);
  for (std::size_t i = 0; i < n; ++i)
    A[smallIndex(i, i)] = 1.0;
}

Real
smallInfinityNorm(const SmallMatrix & A, const std::size_t n)
{
  Real norm = 0.0;

  for (std::size_t i = 0; i < n; ++i)
  {
    Real row_sum = 0.0;
    for (std::size_t j = 0; j < n; ++j)
      row_sum += std::abs(A[smallIndex(i, j)]);

    norm = std::max(norm, row_sum);
  }

  return norm;
}

void
smallMatMul(const SmallMatrix & A,
            const SmallMatrix & B,
            SmallMatrix & C,
            const std::size_t n)
{
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      C[smallIndex(i, j)] = 0.0;

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t k = 0; k < n; ++k)
    {
      const Real aik = A[smallIndex(i, k)];
      if (aik == 0.0)
        continue;

      for (std::size_t j = 0; j < n; ++j)
        C[smallIndex(i, j)] += aik * B[smallIndex(k, j)];
    }
}

SmallExpStats
smallMatrixExponential(const SmallMatrix & A,
                       const std::size_t n,
                       SmallMatrix & result)
{
  if (n == 0 || n > MAX_SMALL_DIM)
    mooseError("MSRPhiFunctions: invalid small matrix dimension.");

  const Real norm = smallInfinityNorm(A, n);

  unsigned int squarings = 0;
  if (norm > 0.5)
    squarings = static_cast<unsigned int>(std::ceil(std::log2(norm / 0.5)));

  const Real scale = std::ldexp(1.0, -static_cast<int>(squarings));

  SmallMatrix B{};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      B[smallIndex(i, j)] = A[smallIndex(i, j)] * scale;

  smallIdentity(result, n);

  SmallMatrix term_a{};
  SmallMatrix term_b{};
  smallIdentity(term_a, n);

  SmallMatrix * term = &term_a;
  SmallMatrix * scratch = &term_b;

  const Real tolerance = 1e-15;
  const unsigned int max_terms = 200;

  bool converged = false;
  unsigned int terms_used = 0;

  for (unsigned int j = 1; j <= max_terms; ++j)
  {
    smallMatMul(*term, B, *scratch, n);

    const Real inv_j = 1.0 / static_cast<Real>(j);

    Real term_norm = 0.0;
    Real result_norm = 0.0;

    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t k = 0; k < n; ++k)
      {
        const auto idx = smallIndex(i, k);
        (*scratch)[idx] *= inv_j;
        result[idx] += (*scratch)[idx];

        term_norm = std::max(term_norm, std::abs((*scratch)[idx]));
        result_norm = std::max(result_norm, std::abs(result[idx]));
      }

    std::swap(term, scratch);
    terms_used = j;

    if (term_norm <= tolerance * std::max(1.0, result_norm))
    {
      converged = true;
      break;
    }
  }

  if (!converged)
    mooseError("MSRPhiFunctions: small scaled matrix exponential Taylor series did not converge.");

  SmallMatrix square_scratch{};

  for (unsigned int i = 0; i < squarings; ++i)
  {
    smallMatMul(result, result, square_scratch, n);

    for (std::size_t r = 0; r < n; ++r)
      for (std::size_t c = 0; c < n; ++c)
        result[smallIndex(r, c)] = square_scratch[smallIndex(r, c)];
  }

  return {terms_used, squarings};
}

std::vector<Real>
phiActionReference(const DenseMatrix & A,
                   const std::vector<Real> & v,
                   const unsigned int k)
{
  const std::size_t n = A.size();
  const std::size_t m = n + k;

  DenseMatrix augmented(m, std::vector<Real>(m, 0.0));

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      augmented[i][j] = A[i][j];

  for (std::size_t i = 0; i < n; ++i)
    augmented[i][n] = v[i];

  for (unsigned int j = 0; j + 1 < k; ++j)
    augmented[n + j][n + j + 1] = 1.0;

  const DenseMatrix exponential = matrixExponential(augmented);

  std::vector<Real> result(n, 0.0);
  const std::size_t column = n + k - 1;

  for (std::size_t i = 0; i < n; ++i)
    result[i] = exponential[i][column];

  return result;
}

std::vector<Real>
phiActionSmall(const DenseMatrix & A,
               const std::vector<Real> & v,
               const unsigned int k)
{
  const std::size_t n = A.size();
  const std::size_t m = n + k;

  SmallMatrix augmented{};

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      augmented[smallIndex(i, j)] = A[i][j];

  for (std::size_t i = 0; i < n; ++i)
    augmented[smallIndex(i, n)] = v[i];

  for (unsigned int j = 0; j + 1 < k; ++j)
    augmented[smallIndex(n + j, n + j + 1)] = 1.0;

  SmallMatrix exponential{};
  smallMatrixExponential(augmented, m, exponential);

  std::vector<Real> result(n, 0.0);
  const std::size_t column = n + k - 1;

  for (std::size_t i = 0; i < n; ++i)
    result[i] = exponential[smallIndex(i, column)];

  return result;
}



void
smallMatVec(const SmallMatrix & A,
            const SmallVector & x,
            SmallVector & y,
            const std::size_t n)
{
  y.fill(0.0);

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      y[i] += A[smallIndex(i, j)] * x[j];
}

void
smallPhiVectorSeries(const SmallMatrix & B,
                     const std::size_t n,
                     const SmallVector & v,
                     const unsigned int k,
                     SmallVector & result)
{
  if (k == 0)
    mooseError("MSRPhiFunctions: phi series index must be at least one.");

  Real factorial = 1.0;
  for (unsigned int j = 2; j <= k; ++j)
    factorial *= static_cast<Real>(j);

  SmallVector term{};
  result.fill(0.0);

  for (std::size_t i = 0; i < n; ++i)
  {
    term[i] = v[i] / factorial;
    result[i] = term[i];
  }

  const Real tolerance = 1e-15;
  const unsigned int max_terms = 200;

  SmallVector scratch{};

  for (unsigned int j = 1; j <= max_terms; ++j)
  {
    smallMatVec(B, term, scratch, n);

    const Real inv = 1.0 / static_cast<Real>(j + k);
    Real term_norm = 0.0;
    Real result_norm = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
      term[i] = scratch[i] * inv;
      result[i] += term[i];

      term_norm = std::max(term_norm, std::abs(term[i]));
      result_norm = std::max(result_norm, std::abs(result[i]));
    }

    if (term_norm <= tolerance * std::max(1.0, result_norm))
      return;
  }

  mooseError("MSRPhiFunctions: direct scaled phi-vector Taylor series did not converge.");
}

std::pair<std::vector<Real>, std::vector<Real>>
phi1Phi3ActionsDirectScaled(const DenseMatrix & A,
                            const std::vector<Real> & v1,
                            const std::vector<Real> & v3)
{
  const std::size_t n = A.size();

  if (n == 0)
    return {{}, {}};

  if (n > MAX_SMALL_DIM)
    mooseError("MSRPhiFunctions: direct scaled phi action exceeds small backend.");

  const Real norm = matrixInfinityNorm(A);

  unsigned int squarings = 0;
  if (norm > 0.5)
    squarings =
        static_cast<unsigned int>(std::ceil(std::log2(norm / 0.5)));

  const Real scale = std::ldexp(1.0, -static_cast<int>(squarings));

  SmallMatrix B{};
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      B[smallIndex(i, j)] = A[i][j] * scale;

  SmallMatrix E{};
  const auto base_exp_stats = smallMatrixExponential(B, n, E);
  if (base_exp_stats.squarings != 0)
    mooseError("MSRPhiFunctions: unexpected secondary scaling in direct phi backend.");

  SmallVector in1{};
  SmallVector in3{};

  for (std::size_t i = 0; i < n; ++i)
  {
    in1[i] = v1[i];
    in3[i] = v3[i];
  }

  SmallVector p1_v1{};
  SmallVector p1_v3{};
  SmallVector p2_v3{};
  SmallVector p3_v3{};

  smallPhiVectorSeries(B, n, in1, 1, p1_v1);
  smallPhiVectorSeries(B, n, in3, 1, p1_v3);
  smallPhiVectorSeries(B, n, in3, 2, p2_v3);
  smallPhiVectorSeries(B, n, in3, 3, p3_v3);

  SmallVector Ep1_v1{};
  SmallVector Ep1_v3{};
  SmallVector Ep2_v3{};
  SmallVector Ep3_v3{};

  SmallVector new_p1_v1{};
  SmallVector new_p1_v3{};
  SmallVector new_p2_v3{};
  SmallVector new_p3_v3{};

  SmallMatrix E_squared{};

  for (unsigned int level = 0; level < squarings; ++level)
  {
    smallMatVec(E, p1_v1, Ep1_v1, n);
    smallMatVec(E, p1_v3, Ep1_v3, n);
    smallMatVec(E, p2_v3, Ep2_v3, n);
    smallMatVec(E, p3_v3, Ep3_v3, n);

    for (std::size_t i = 0; i < n; ++i)
    {
      new_p1_v1[i] = 0.5 * (p1_v1[i] + Ep1_v1[i]);
      new_p1_v3[i] = 0.5 * (p1_v3[i] + Ep1_v3[i]);

      new_p2_v3[i] =
          0.25 * (p1_v3[i] + p2_v3[i] + Ep2_v3[i]);

      new_p3_v3[i] =
          0.125 *
          (p2_v3[i] + 0.5 * p1_v3[i] + p3_v3[i] + Ep3_v3[i]);
    }

    p1_v1 = new_p1_v1;
    p1_v3 = new_p1_v3;
    p2_v3 = new_p2_v3;
    p3_v3 = new_p3_v3;

    smallMatMul(E, E, E_squared, n);
    E = E_squared;
  }

  std::vector<Real> out1(n, 0.0);
  std::vector<Real> out3(n, 0.0);

  for (std::size_t i = 0; i < n; ++i)
  {
    out1[i] = p1_v1[i];
    out3[i] = p3_v3[i];
  }

  return {std::move(out1), std::move(out3)};
}

} // namespace

std::vector<Real>
MSRPhiFunctions::matVec(const DenseMatrix & A,
                        const std::vector<Real> & v)
{
  const std::size_t n = A.size();

  if (v.size() != n)
    mooseError("MSRPhiFunctions: matrix/vector size mismatch.");

  std::vector<Real> result(n, 0.0);

  for (std::size_t i = 0; i < n; ++i)
  {
    if (A[i].size() != n)
      mooseError("MSRPhiFunctions: matrix must be square.");

    for (std::size_t j = 0; j < n; ++j)
      result[i] += A[i][j] * v[j];
  }

  return result;
}

std::vector<Real>
MSRPhiFunctions::phiAction(const DenseMatrix & A,
                           const std::vector<Real> & v,
                           unsigned int k)
{
  if (k == 0)
    mooseError("MSRPhiFunctions: phi index must be at least one.");

  const std::size_t n = A.size();

  if (v.size() != n)
    mooseError("MSRPhiFunctions: matrix/vector size mismatch.");

  for (const auto & row : A)
    if (row.size() != n)
      mooseError("MSRPhiFunctions: matrix must be square.");

  const std::size_t m = n + k;

  if (m <= MAX_SMALL_DIM)
  {
    auto result = phiActionSmall(A, v, k);

    static std::atomic<unsigned int> verification_calls{0};
    const auto ticket = verification_calls.fetch_add(1, std::memory_order_relaxed);

    if (ticket < 24)
    {
      const auto reference = phiActionReference(A, v, k);

      for (std::size_t i = 0; i < n; ++i)
      {
        const Real abs_diff = std::abs(result[i] - reference[i]);
        const Real scale = std::max({1.0, std::abs(result[i]), std::abs(reference[i])});

        if (abs_diff > 5e-11 * scale)
          mooseError("MSRPhiFunctions: optimized small-matrix phi action failed reference check.");
      }
    }

    return result;
  }

  return phiActionReference(A, v, k);
}


std::pair<std::vector<Real>, std::vector<Real>>
MSRPhiFunctions::phi1Phi3Actions(const DenseMatrix & A,
                                 const std::vector<Real> & v1,
                                 const std::vector<Real> & v3)
{
  const std::size_t n = A.size();

  if (v1.size() != n || v3.size() != n)
    mooseError("MSRPhiFunctions: combined phi1/phi3 vector size mismatch.");

  for (const auto & row : A)
    if (row.size() != n)
      mooseError("MSRPhiFunctions: matrix must be square.");

  std::pair<std::vector<Real>, std::vector<Real>> result;

  if (n <= MAX_SMALL_DIM)
    result = phi1Phi3ActionsDirectScaled(A, v1, v3);
  else
    result = {phiActionReference(A, v1, 1), phiActionReference(A, v3, 3)};

  static std::atomic<unsigned int> combined_verification_calls{0};
  const auto ticket =
      combined_verification_calls.fetch_add(1, std::memory_order_relaxed);

  if (ticket < 24 || ticket % 10000 == 0)
  {
    const auto ref1 = phiActionReference(A, v1, 1);
    const auto ref3 = phiActionReference(A, v3, 3);

    for (std::size_t i = 0; i < n; ++i)
    {
      const Real diff1 = std::abs(result.first[i] - ref1[i]);
      const Real scale1 =
          std::max({1.0, std::abs(result.first[i]), std::abs(ref1[i])});

      const Real diff3 = std::abs(result.second[i] - ref3[i]);
      const Real scale3 =
          std::max({1.0, std::abs(result.second[i]), std::abs(ref3[i])});

      if (diff1 > 5e-11 * scale1 || diff3 > 5e-11 * scale3)
        mooseError("MSRPhiFunctions: combined phi1/phi3 action failed reference check.");
    }
  }

  return result;
}

std::vector<Real>
MSRPhiFunctions::phi1Action(const DenseMatrix & A,
                            const std::vector<Real> & v)
{
  return phiAction(A, v, 1);
}

std::vector<Real>
MSRPhiFunctions::phi3Action(const DenseMatrix & A,
                            const std::vector<Real> & v)
{
  return phiAction(A, v, 3);
}

} // namespace MSR
