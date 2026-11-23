"""Exploratory per-gene NB mixture with a library-size offset."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import digamma, expit, gammaln, logsumexp


def fit_nb_split(counts, libraries):
    """Inputs contain only positive-depth cells. Return summary and high-state mask."""
    if len(counts) < 4 or np.count_nonzero(counts) == 0:
        return {"status": "insufficient_signal"}, None
    y = np.asarray(counts, dtype=float)[:, None]
    exposure = np.asarray(libraries, dtype=float)[:, None] / 10000
    rates = y[:, 0] / exposure[:, 0]
    mean = y.sum() / exposure.sum()
    log_factorial = gammaln(y + 1)

    def fit(means, fraction):
        k = len(means)

        def objective(params, posterior=False):
            mu = np.exp(params[:k]) * exposure
            r = np.exp(params[-1])
            logits = np.r_[0., params[k:-1]]
            log_weights = logits - logsumexp(logits)
            log_ratio = np.log1p(mu / r)
            joint = (gammaln(y + r) - gammaln(r) - log_factorial
                     + y * np.log(mu / r) - (y + r) * log_ratio + log_weights)
            total = logsumexp(joint, axis=1, keepdims=True)
            weights = np.exp(joint - total)
            if posterior:
                return weights
            grad_mean = (weights * (y - mu) * r / (mu + r)).sum(axis=0)
            grad_logits = weights.sum(axis=0) - len(y) * np.exp(log_weights)
            grad_r = (weights * r * (digamma(y + r) - digamma(r) - log_ratio
                                     + (mu - y) / (mu + r))).sum()
            return -total.sum(), -np.r_[grad_mean, grad_logits[1:], grad_r]

        initial = np.r_[np.log(np.maximum(means, 1e-6)),
                        [np.log(fraction / (1 - fraction))] * (k - 1), 0.]
        bounds = [(-16, 14)] * k + [(-20, 20)] * (k - 1) + [(-12, 9)]
        result = minimize(objective, initial, jac=True, method="L-BFGS-B", bounds=bounds,
                          options={"ftol": 1e-12})
        return result, objective

    single, _ = fit([mean], 0.5)
    fits = []
    for q in (0.5, 0.95, 0.998):
        cut = np.quantile(rates, q)
        upper = rates > cut
        if not upper.any():
            upper = rates == rates.max()
        if upper.all():
            continue
        means = [y[mask].sum() / exposure[mask].sum() for mask in (~upper, upper)]
        result, objective = fit(means, upper.mean())
        if result.success and np.isfinite(result.fun):
            fits.append((result, objective))
    if not fits or not single.success:
        return {"status": "no_reliable_fit"}, None
    best, objective = min(fits, key=lambda pair: pair[0].fun)
    delta = 2 * (single.fun - best.fun - np.log(len(y)))
    means = np.exp(best.x[:2])
    high = int(np.argmax(means))
    summary = {"status": "two_components" if delta > 0 else "one_component",
               "bic_one_minus_two": float(delta), "means_per_10k": sorted(means.tolist()),
               "high_component_weight": float(expit(best.x[2] if high else -best.x[2]))}
    if delta <= 0:
        return summary, None
    return summary, objective(best.x, posterior=True)[:, high] > 0.5
