# Polypharmacy Side-Effect Burden

Predicting side effects of drugs taken together, and weighting those predictions by how often people actually take them together.

## Language

**Drug pair**:
Exactly two distinct drugs taken together; the unit every prediction, prevalence and burden figure is about.
_Avoid_: Combination, combo, interaction

**Drug set**:
Three or more drugs taken together, such as one survey participant's full medication list; scored only by combining its drug pairs' scores.
_Avoid_: Combination, regimen, cocktail

**Labelled pair**:
A drug pair with at least one side effect in the Decagon data.
_Avoid_: Positive pair, interacting pair, known interaction

**Unlabelled pair**:
A drug pair with no side effect in the Decagon data; this means "not reported", not "safe".
_Avoid_: Negative pair, safe pair, non-interacting pair
