# Agent Handoff Notes

## User preferences for `qe_convergence_tests.ipynb`

- Keep the notebook very close to the original exercise text and workflow.
- Use simple language with minimal notes.
- Leave space for student action instead of over-explaining.
- Build QE inputs explicitly with `pw_input.py` namelist/card objects.
- Prefer pythonic sweep interfaces and avoid explicit `for` loops where practical.
- Keep plotting/helper material modular or outside the main exercise notebook when possible.
- Include optional force/stress convergence.
- For force convergence, explicitly include displacement of one atom first.

## Session context

- Intro rationale was added at the notebook start about why convergence tests are needed.
- User asked to stop and continue with a different agent.
