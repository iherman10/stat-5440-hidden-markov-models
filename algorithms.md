# Adapting the Lecture HMM to the Multi-Game Baseball Example

## Big Picture

The core algorithms (FFBS, Gibbs sampler) are unchanged in their logic. The adaptation from the lecture's single-sequence, 2-state HMM to the baseball setting boils down to two things:

1. **Wrap FFBS in a loop over games**, treating each game as an independent sequence.
2. **Aggregate counts carefully** across games when updating shared parameters.

This document covers **four permutations** of the baseball model:

|                                   | Raw Probabilities | Log-Space |
| --------------------------------- | ----------------- | --------- |
| **2-State** (Baseline, Hot)       | Section A         | Section B |
| **3-State** (Cold, Baseline, Hot) | Section C         | Section D |

All four share the same game-loop structure and Gibbs sampler skeleton. What changes across them is the dimension of the state space (affecting priors and label switching) and whether the forward pass uses raw or log-scale arithmetic (affecting only the alpha table computation inside FFBS).

---

## Shared Structure: The Game Loop

Regardless of the number of states or the numerical scale, the multi-game structure works the same way.

We have $G$ games (starts), each of length $T_g$ (number of four-seam fastballs in that game). Because hidden states **reset at the start of each game**, each game is an independent sequence. Per Gibbs iteration, FFBS runs $G$ times. Each run:

- **Initializes** the alpha table from $\pi$ (the shared initial state distribution), not from the previous game's final state.
- **Fills** the alpha table forward using $\nu$ and the emission densities.
- **Samples** the state sequence backward.

There is **no transition** between the last pitch of game $g$ and the first pitch of game $g+1$.

---

## Section A: 2-State Model, Raw Probabilities

This is the closest to the lecture code. You are essentially wrapping the lecture's single-sequence logic in a game loop and swapping the priors.

### Priors

```r
K <- 2

prior_nu <- rbind(
    c(8, 1),  # Baseline row
    c(1, 8)   # Hot row
)

prior_pi <- c(6, 2)  # Baseline-heavy start
```

### FFBS (Block 1)

```r
for (g in 1:G) {

    # ---- FORWARD PASS ----
    for (i in 1:K) {
        alpha[1, i] <- curpi[i] * dnorm(Y[[g]][1], curmu[i], sigma)
    }

    for (t in 2:T[g]) {
        for (i in 1:K) {
            alpha[t, i] <- 0
            for (j in 1:K) {
                alpha[t, i] <- alpha[t, i] +
                    alpha[t-1, j] * curnu[j, i] * dnorm(Y[[g]][t], curmu[i], sigma)
            }
        }
    }

    # ---- BACKWARD PASS ----
    prob_vec <- alpha[T[g], ] / sum(alpha[T[g], ])
    X[[g]][T[g]] <- sample(1:K, size = 1, prob = prob_vec)

    for (t in (T[g] - 1):1) {
        for (i in 1:K) {
            prob_vec[i] <- alpha[t, i] * curnu[i, X[[g]][t+1]]
        }
        prob_vec <- prob_vec / sum(prob_vec)
        X[[g]][t] <- sample(1:K, size = 1, prob = prob_vec)
    }
}
```

### Transition Matrix Update (Block 2)

```r
Nmat <- matrix(0, nrow = K, ncol = K)
for (g in 1:G) {
    for (t in 2:T[g]) {
        j <- X[[g]][t-1]
        k <- X[[g]][t]
        Nmat[j, k] <- Nmat[j, k] + 1
    }
}

for (j in 1:K) {
    for (k in 1:K) {
        curnu[j, k] <- rgamma(1, shape = Nmat[j, k] + prior_nu[j, k], rate = 1)
    }
    curnu[j, ] <- curnu[j, ] / sum(curnu[j, ])
}
```

### Initial Distribution Update (Block 3)

```r
init_count <- rep(0, K)
for (g in 1:G) {
    init_count[X[[g]][1]] <- init_count[X[[g]][1]] + 1
}

for (j in 1:K) {
    curpi[j] <- rgamma(1, shape = init_count[j] + prior_pi[j], rate = 1)
}
curpi <- curpi / sum(curpi)
```

### Emission Means Update (Block 4)

```r
for (j in 1:K) {
    # Pool all pitches assigned to state j across all games
    obs_j <- unlist(lapply(1:G, function(g) Y[[g]][X[[g]] == j]))
    n_j <- length(obs_j)

    if (n_j > 0) {
        Y_bar_j <- mean(obs_j)
        curmu[j] <- rnorm(1, mean = Y_bar_j, sd = sqrt(sigma^2 / n_j))
    }
}
```

### Label Switching

Same as lecture: if $\mu_1 > \mu_2$, swap $\mu$, swap $\pi$, permute rows and columns of $\nu$, and relabel all state sequences.

```r
if (curmu[1] > curmu[2]) {
    curmu <- rev(curmu)
    curpi <- rev(curpi)
    curnu <- curnu[2:1, 2:1]
    for (g in 1:G) {
        X[[g]] <- ifelse(X[[g]] == 1, 2, 1)
    }
}
```

---

## Section B: 2-State Model, Log-Space

Identical to Section A in all Gibbs blocks (2–4) and label switching. The only change is inside FFBS, where the alpha table is computed in log-space using the **log-sum-exp trick**.

### Why log-space?

Multiplying many probabilities < 1 can underflow to zero. Log-space replaces products with sums, keeping values numerically stable. The log-sum-exp trick handles the summation step: $\log(\sum_j e^{a_j}) = a_{\max} + \log(\sum_j e^{a_j - a_{\max}})$.

### FFBS in Log-Space (Block 1)

```r
for (g in 1:G) {

    # ---- FORWARD PASS (log-space) ----
    for (i in 1:K) {
        log_alpha[1, i] <- log(curpi[i]) + dnorm(Y[[g]][1], curmu[i], sigma, log = TRUE)
    }

    for (t in 2:T[g]) {
        for (i in 1:K) {
            log_candidates <- rep(0, K)
            for (j in 1:K) {
                log_candidates[j] <- log_alpha[t-1, j] +
                    log(curnu[j, i]) +
                    dnorm(Y[[g]][t], curmu[i], sigma, log = TRUE)
            }
            # Log-sum-exp trick
            max_val <- max(log_candidates)
            log_alpha[t, i] <- max_val + log(sum(exp(log_candidates - max_val)))
        }
    }

    # ---- BACKWARD PASS (pop out of log-space only to sample) ----
    log_probs <- log_alpha[T[g], ]
    max_val <- max(log_probs)
    prob_vec <- exp(log_probs - max_val)
    prob_vec <- prob_vec / sum(prob_vec)
    X[[g]][T[g]] <- sample(1:K, size = 1, prob = prob_vec)

    for (t in (T[g] - 1):1) {
        log_prob_vec <- rep(0, K)
        for (i in 1:K) {
            log_prob_vec[i] <- log_alpha[t, i] + log(curnu[i, X[[g]][t+1]])
        }
        max_val <- max(log_prob_vec)
        prob_vec <- exp(log_prob_vec - max_val)
        prob_vec <- prob_vec / sum(prob_vec)
        X[[g]][t] <- sample(1:K, size = 1, prob = prob_vec)
    }
}
```

Blocks 2–4 and label switching are **identical to Section A**. Log-space only affects the alpha table inside FFBS.

---

## Section C: 3-State Model, Raw Probabilities

The 3-state model adds the "Cold/Fatigued" state. The FFBS and Gibbs block logic is the same as Section A, with all loops running over `1:3` instead of `1:2` and updated priors.

### Priors

```r
K <- 3

prior_nu <- rbind(
    c(8, 1, 1),  # Cold row
    c(1, 8, 1),  # Baseline row
    c(1, 1, 8)   # Hot row
)

prior_pi <- c(2, 6, 2)  # Baseline-heavy start
```

### FFBS (Block 1)

Identical structure to Section A, but all inner loops iterate over `1:3` instead of `1:2`. No other changes.

### Blocks 2–4

Identical structure to Section A with `K <- 3`. The transition count matrix is now $3 \times 3$, the initial count vector has 3 entries, and emission means are estimated for 3 states.

### Label Switching (3 States)

This is the key implementation difference from the 2-state case. With $3! = 6$ possible permutations, the simple swap no longer works. Instead, find the sort permutation and apply it everywhere:

```r
perm <- order(curmu)  # e.g., c(2, 3, 1) if mu[2] < mu[3] < mu[1]

# Apply permutation
curmu <- curmu[perm]
curpi <- curpi[perm]
curnu <- curnu[perm, perm]  # permute BOTH rows and columns

for (g in 1:G) {
    # Build inverse mapping: inv_perm[old_label] = new_label
    inv_perm <- rep(0, K)
    inv_perm[perm] <- 1:K
    X[[g]] <- inv_perm[X[[g]]]
}
```

The `order()` call returns the permutation that sorts $\mu$ into $\mu_1 < \mu_2 < \mu_3$ (Cold < Baseline < Hot). Applying it to rows **and** columns of $\nu$ is critical — rows handle "from" state relabeling, columns handle "to" state relabeling.

---

## Section D: 3-State Model, Log-Space

Combines the 3-state priors and label switching from Section C with the log-space FFBS from Section B.

### FFBS in Log-Space (Block 1)

```r
for (g in 1:G) {

    # ---- FORWARD PASS (log-space) ----
    for (i in 1:K) {
        log_alpha[1, i] <- log(curpi[i]) + dnorm(Y[[g]][1], curmu[i], sigma, log = TRUE)
    }

    for (t in 2:T[g]) {
        for (i in 1:K) {
            log_candidates <- rep(0, K)
            for (j in 1:K) {
                log_candidates[j] <- log_alpha[t-1, j] +
                    log(curnu[j, i]) +
                    dnorm(Y[[g]][t], curmu[i], sigma, log = TRUE)
            }
            max_val <- max(log_candidates)
            log_alpha[t, i] <- max_val + log(sum(exp(log_candidates - max_val)))
        }
    }

    # ---- BACKWARD PASS ----
    log_probs <- log_alpha[T[g], ]
    max_val <- max(log_probs)
    prob_vec <- exp(log_probs - max_val)
    prob_vec <- prob_vec / sum(prob_vec)
    X[[g]][T[g]] <- sample(1:K, size = 1, prob = prob_vec)

    for (t in (T[g] - 1):1) {
        log_prob_vec <- rep(0, K)
        for (i in 1:K) {
            log_prob_vec[i] <- log_alpha[t, i] + log(curnu[i, X[[g]][t+1]])
        }
        max_val <- max(log_prob_vec)
        prob_vec <- exp(log_prob_vec - max_val)
        prob_vec <- prob_vec / sum(prob_vec)
        X[[g]][t] <- sample(1:K, size = 1, prob = prob_vec)
    }
}
```

### Blocks 2–4

Identical to Section C (3-state, raw space). Log-space only affects the alpha table inside FFBS.

### Label Switching

Identical to Section C (use `order()` + full permutation logic).

---

## Implementation Gotchas

### Numerical Underflow (Sections A and C)

Raw-space implementations can underflow to zero for longer sequences or tighter emission distributions. Symptoms: `NaN`s or all-zero alpha rows. Fix: switch to log-space (Sections B/D).

With ~30–60 pitches per game and 2 states, you'll likely be fine in raw space. With 3 states or if you ever extend to longer sequences, log-space is the safer default.

### Label Switching Complexity

| States | Permutations | Fix                                                                                               |
| ------ | ------------ | ------------------------------------------------------------------------------------------------- |
| 2      | 2            | Simple swap: if $\mu_1 > \mu_2$, swap everything                                                  |
| 3      | 6            | Sort $\mu$ via `order()`, apply permutation to $\mu$, $\pi$, rows+cols of $\nu$, and state labels |

Getting the 3-state permutation wrong (e.g., permuting only rows of $\nu$ but not columns) will silently corrupt your posterior. Test by checking that $\mu_1 < \mu_2 < \mu_3$ holds for every post-processed sample.

### Edge Case: Empty States

If a state has zero observations assigned to it in a Gibbs iteration (possible early on or with 3 states), the emission mean update will fail (division by zero). Handle this by keeping the previous value:

```r
if (n_j > 0) {
    Y_bar_j <- mean(obs_j)
    curmu[j] <- rnorm(1, mean = Y_bar_j, sd = sqrt(sigma^2 / n_j))
}
# else: curmu[j] stays unchanged from previous iteration
```

---

## Quick Reference: What Changes Across Permutations

| Component                  | 2-State Raw (A)       | 2-State Log (B)           | 3-State Raw (C)                    | 3-State Log (D)           |
| -------------------------- | --------------------- | ------------------------- | ---------------------------------- | ------------------------- |
| **FFBS alpha table**       | Raw products          | Log sums + log-sum-exp    | Raw products                       | Log sums + log-sum-exp    |
| **FFBS backward sampling** | Normalize raw alpha   | exp(log - max), normalize | Normalize raw alpha                | exp(log - max), normalize |
| **Transition prior**       | Dirichlet(8,1), (1,8) | Same as A                 | Dirichlet(8,1,1), (1,8,1), (1,1,8) | Same as C                 |
| **$\pi$ prior**            | Dirichlet(6, 2)       | Same as A                 | Dirichlet(2, 6, 2)                 | Same as C                 |
| **Label switching**        | Simple swap           | Same as A                 | `order()` + full permutation       | Same as C                 |
| **Gibbs Blocks 2–4**       | $K = 2$               | Same as A                 | $K = 3$                            | Same as C                 |
| **Underflow risk**         | Low                   | None                      | Moderate                           | None                      |
