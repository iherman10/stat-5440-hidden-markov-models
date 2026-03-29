## STAT 5440 Final Project

## Applied Bayesian Mini-Lecture Project

## Overview

In groups, you will design and deliver a $\mathbf{4 0}$-minute mini-lecture on an applied Bayesian topic. Your goal is to:

1. Motivate the method
2. Apply it to a real dataset
3. Teach the class how to implement it
4. Provide a clean, reproducible tutorial in a markdown format.

This project simulates what it means to teach and communicate applied Bayesian methods clearly and responsibly.

## Project Components

Each group must produce:

## 1 In-Class Mini-Lecture ( $\mathbf{4 0}$ minutes $\mathbf{+} \mathbf{5 ~ m i n}$ Q\&A)

The lecture must include:

## A. Motivation

- What problem does this method solve?
- Why can't standard regression handle it?
- Where does it show up in practice?


## B. Model Specification

- Likelihood
- Priors
- Hierarchical structure (if relevant)
- Key modeling assumptions


## C. Applied Example

- Real dataset (not simulated-only)
- Clear explanation of data context
- Model fit using stan / rstanarm / cmdstandr / brms / PyStan / CmdStanPy
- Diagnostics
- Posterior summaries
- Posterior predictive checks


## D. Interpretation

- What did we learn?
- What are limitations?
- When should we use this model?

Every student must speak during the lecture. Dominance of speaking time of a subset of students in your group will negatively impact your grade.

## 2 Tutorial Document (Reproducible)

Each group must submit a tutorial document that:

- Explains the method in written form
- Includes clean, well-commented code
- Runs start-to-finish without error
- Produces the figures shown in class

Acceptable formats:

- RMarkdown (.Rmd)
- Quarto (.qmd)
- Jupyter Notebook (.ipynb)

The script that creates the tutorial (.Rmd/.qmd/.ipynb) and the corresponding .html document will need to be uploaded to the corresponding class github page for projects by the due date in order to be considered completed.

The tutorial must:

- Explain the motivation behind the problem
- Give and explain the model specification
- Load data
- Fit model
- Run diagnostics
- Generate relevant plots
- Include interpretation

Every student must contribute to the tutorial. Every student must respond to a questionnaire about how much each student contributed to the tutorial/project work. A hierarchical model will be used to assign credit based on such responses (to be specified later).

## Grading Rubric (100 points)

## Presentation (50 points)

## Category

Conceptual clarity
Model explanation
Applied example quality

## Points

10
10
10
10
5
5

## Tutorial Document (50 points)

## Category

Compelling Motivation

## Points

5

Explanation of model specification 5
Correct implementation 5
Reproducibility 10
Code clarity \& commenting 10
Interpretation of results 10
Professional formatting5

## Group Structure \& Scheduling

There are:

- 26 students
- 40 min presentation +5 min Q\&A = 45 min per group
- 90 min per class
- Up to 3 class sessions available


## Plan

- 6 groups of 4-5 students
- 2 groups per class
- 3 class sessions total

Each class: Two 45 minute presentations (90 minutes of class total)

## Topic List

Students will rank preferences (1-6). You assign topics using ranked matching.

## Core Topics

1. Gaussian Processes
2. Nonparametric Bayesian Methods/Dirichlet Processes
3. Spatial Models
4. Dynamic Linear Models / State-space Models / Time Series
5. Survival Models
6. Censored \& Missing Data
7. Newer Approximation / Sampling Methods (Pathfinder, Laplace, INLA, etc.)
8. Hidden Markov Models (HMM)
9. Hurdle Models (and other zero-inflated models not covered in class)
10. Causal Inference with Bayesian Methods
11. Bayesian Networks
12. Change Point Models
13. Measurement Error Models

## Topic Assignment Mechanism

To fairly separate groups based on ranked preferences:

## Collect Preferences

Each student submits:

- Ranked top 6 topics
- Optional: preferred collaborator (max 1 request)


## Form Groups \& Assign Topics

1. I will try to pair people who requested each other first
2. Assign topics in order of random group priority.
3. Give each group their highest-ranked available topic.
4. Remove that topic from pool.
5. Continue until all assigned.

## Timeline

- Preferences submitted by Friday March $6{ }^{\text {th }}$, 2026
- Form: https://forms.gle/ngLfgoLq7QNB7bHg6
- Groups assigned shortly thereafter
- Presentations start Tuesday April 21, 2026
- Tutorials uploaded and submitted by Wednesday April $29{ }^{\text {th }}$, by 11:59 PM


## Requirements \& Constraints

In order to have a great project you should

## Include:

- Real dataset (not toy-only)
- Posterior predictive check
- MCMC diagnostics
- At least one meaningful visualization
- Interpretation in context
- How would results change under a different assumption?
- What is one modeling decision that materially affects inference?


## And Cannot:

- Just restate textbook example
- Use pre-built black-box without explanation
- Skip model assumptions


## What this project is testing

- Can you translate theory → practice?
- Can you explain modeling assumptions?
- Can you produce reproducible research?
- Do you understand posterior interpretation?

