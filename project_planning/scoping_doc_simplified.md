# Simplified Scoping Document: Bayesian HMM for Pitcher Performance

## Project Overview

We model latent performance modes of a single MLB starting pitcher using a Bayesian Hidden Markov Model. The idea is that a pitcher's "true condition" during a game is not directly observable — we only see noisy pitch-level measurements. A pitcher might be locked in for 30 pitches, lose velocity for 15 pitches, then stabilize. An HMM captures this kind of regime-switching behavior naturally.

We observe pitch-level velocity data for four-seam fastballs across one season of starts. The model infers when the pitcher is likely in a "hot," "baseline," or "cold/fatigued" state, and how those states evolve within and across games.

## Key Simplifications (vs. Original Scope)

- **One pitcher, one season** (no hierarchical/multi-pitcher structure)
- **One emission metric: velocity** (univariate Gaussian, not multivariate)
- **Standard Dirichlet priors** on transition rows (no sticky Dirichlet)
- **No covariates** (e.g., pitch count, inning) — pure HMM to start

## Data

- **Source:** MLB Statcast via Baseball Savant or the `pybaseball` / `baseballr` packages
- **Granularity:** Pitch-level data for a single pitcher across one full season
- **Filter:** Four-seam fastballs only
  - Different pitch types (slider, changeup, curveball) have inherently different velocity and movement profiles. If we mix pitch types, the HMM would likely discover "pitch type" as the hidden state rather than performance mode. Restricting to four-seam fastballs keeps the hidden state interpretation clean.
- **Key field:** `release_speed` (pitch velocity in mph)
- **Pitcher selection:** TBD — to be decided by group. Should be a starter with high four-seam fastball usage and a full season of starts (~30+ starts). This gives us roughly 1,000–2,000 four-seam fastballs to work with.
- **Note on emission metric:** We are starting with velocity as the single emission metric because it is the most interpretable and well-studied in this context. This can be swapped for another metric (e.g., spin rate, release extension, or induced vertical break) if the group prefers, without changing the model structure.

## Model Specification

### Structure

- Let $g = 1, \dots, G$ index games (starts) in the season
- Let $t = 1, \dots, T_g$ index four-seam fastball pitches within game $g$
- Each game is treated as an **independent sequence** — the hidden state resets at the start of each game (a pitcher's fatigue state in one game tells us nothing about the next game, which is typically 5 days later)

### Hidden States

$$z_{g,t} \in \{1, 2, 3\}$$

- **State 1: Cold/Fatigued** — lower velocity
- **State 2: Baseline** — typical velocity for this pitcher
- **State 3: Hot** — elevated velocity

### Emission Model

Within each state, observed velocity is modeled as a univariate Gaussian:

$$y_{g,t} \mid z_{g,t} = k \sim \text{Normal}(\mu_k, \sigma_k^2)$$

where $\mu_k$ is the state-specific mean velocity and $\sigma_k^2$ is the state-specific variance.

### Transition Model

The pitcher's hidden state at each pitch depends on the state at the previous pitch. The relationship is governed by a **transition matrix** $A$, which is a $3 \times 3$ table of probabilities:

| | Next: Cold | Next: Baseline | Next: Hot |
|---|---|---|---|
| **Current: Cold** | high | low | low |
| **Current: Baseline** | low | high | low |
| **Current: Hot** | low | low | high |

Each row sums to 1. To determine the next state, we look up the current state's row and draw from those probabilities. For example, if the pitcher is currently in the Baseline state, we look at row 2 and draw the next state from those three probabilities. In formal notation:

$$z_{g,t} \mid z_{g,t-1} \sim \text{Categorical}(A[z_{g,t-1}, \cdot])$$

This just means: "pick the next state randomly according to the probabilities in the current state's row."

### Initial State Distribution

At the start of each game:

$$z_{g,1} \sim \text{Categorical}(\pi)$$

where $\pi$ is the initial state probability vector.

## Priors

### Emission Means

Weakly informative normal priors on each state's mean velocity:

$$\mu_k \sim \text{Normal}(\bar{y}, \, 5^2)$$

where $\bar{y}$ is the pitcher's overall average four-seam velocity (computed from data). The variance of 25 is deliberately wide — it says we expect state means to be somewhere in the neighborhood of the pitcher's average, but we're not being very prescriptive.

### Emission Variances

We need a prior on $\sigma_k$ (the state-specific standard deviation). Two common options:

1. **Exponential prior (simplest):**
   $$\sigma_k \sim \text{Exponential}(1)$$
   Simple, one parameter, weakly informative. Puts most mass on small-to-moderate standard deviations. This is the recommended default for simplicity.

2. **Half-Cauchy prior:**
   $$\sigma_k \sim \text{Half-Cauchy}(0, 2.5)$$
   Slightly heavier tails than the exponential, which makes it more forgiving if the data has more spread than expected. A standard choice in Bayesian modeling but marginally more complex.

**Recommendation:** Start with the Exponential(1) prior. It's the simplest and works well for this kind of problem.

### Transition Matrix Rows

Each row of the transition matrix gets a Dirichlet prior. The Dirichlet parameters control what transition probabilities we expect before seeing data. Larger values mean "we expect more probability here."

We want each state to favor staying in itself (persistence), so we put a large value on the diagonal entry and small values elsewhere:

| Row (current state) | Dirichlet prior | Interpretation |
|---|---|---|
| **Cold** | Dirichlet(8, 1, 1) | Expect ~80% chance of staying cold |
| **Baseline** | Dirichlet(1, 8, 1) | Expect ~80% chance of staying baseline |
| **Hot** | Dirichlet(1, 1, 8) | Expect ~80% chance of staying hot |

The "8" on the diagonal encodes our belief that states persist. The "1"s on the off-diagonal say transitions are possible but less likely. These are soft priors — the data can easily override them. The exact values (e.g., 8 vs. 6 vs. 10) can be tuned.

### Initial State Distribution

$$\pi \sim \text{Dirichlet}(2, 6, 2)$$

| State | Dirichlet parameter | Prior expected probability |
|---|---|---|
| Cold | 2 | ~20% |
| Baseline | 6 | ~60% |
| Hot | 2 | ~20% |

This says we expect most games to begin near the baseline state, with some probability of starting hot or cold. This is a soft prior — the data can easily move it.

## Core Assumptions

1. **First-order Markov:** The current hidden state depends only on the previous state, not the full history.
2. **Conditional independence:** Given the hidden state, pitch velocities are independent across time steps.
3. **Stationary transitions:** The same transition matrix applies throughout the season (no game-to-game drift in transition dynamics).
4. **Game independence:** Hidden states reset at the start of each game.

## Possible Extensions (Not in Scope for v1)

These are things we could add later if time permits or the base model works well:

- Add a second emission metric (e.g., spin rate) and move to multivariate Gaussian emissions
- Add covariates to the transition matrix (e.g., pitch count, inning) to model fatigue effects
- Expand to multiple pitchers with hierarchical priors on emission parameters
- Use sticky Dirichlet priors instead of standard Dirichlet for stronger persistence
- Compare 2-state vs. 3-state models using model comparison criteria (e.g., WAIC, LOO-CV)
