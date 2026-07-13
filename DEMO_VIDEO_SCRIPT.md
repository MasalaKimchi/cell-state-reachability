# 3-Minute Demo Video — Shooting Script

**Project:** Cell-State Reachability

**Event:** Built with Claude — Life Sciences Hackathon (Research / Lab Track)
**Runtime target:** 2:45–3:00 · **Voiceover:** ~360 words
**Published video:** https://www.youtube.com/watch?v=GJbxLxYUBMo

## The one-sentence hook

> **A GPS does not just suggest a route—it tells you when the road you need is not on the
> map. We built that second answer for cell engineering.**

## Beat sheet

| # | Time | On screen | Voiceover |
|---|---:|---|---|
| 1 | 0:00–0:20 | `app/explorers/pharma_funnel.html`, then move quickly to the “target selection” box | “Drug discovery pays for wrong turns late. This tool does not claim to predict clinical success. It asks a smaller question early: with the intervention effects we have actually measured, can we even point the cell in the direction we want?” |
| 2 | 0:20–0:42 | `fig_central_illustration.png`, panel A | “Start with a real genome-scale CRISPRi screen in primary human T cells: thirty-three thousand nine hundred eighty-three perturbation-by-condition profiles across more than ten thousand readout genes. Each perturbation leaves a measured fingerprint—some genes rise, others fall.” |
| 3 | 0:42–1:08 | `app/explorers/reachability_explorer.html`, geometry panel | “We may combine those fingerprints, but we cannot reverse one and pretend we performed a negative knockdown. Under a first-order additivity model, all non-negative combinations form a cone. The target state is an arrow. The question is: how closely can this measured cone point along it?” |
| 4 | 1:08–1:38 | Explorer decomposition, then certificate panel | “For the Th2-to-Th1 direction, the answer is: partly. The held-out cosine is point four four eight, higher than all sixty shuffled targets. Thirty-nine percent of target energy is captured by measured knockdown directions. Twenty-five percent is captured only by a sign-flipped activation proxy, and thirty-five percent by neither. The residual proves the full target arrow sits outside the measured knockdown cone.” |
| 5 | 1:38–2:02 | Certificate genes, then `pharma_triage.html` | “Its largest positive coordinates rank what the closest knockdown mixture still under-delivers. Those genes are not proven activation targets—they are a short list for CRISPRa or de-repression tests. Separately, the top greedy knockdown panels collapse to one hundred two unique candidates; forty-five lack a conventional drug handle and ten have a clinical-grade drug in the saved snapshot.” |
| 6 | 2:02–2:24 | `pharma_capability.html` | “We structured a survey of ninety-one prior methods. Some use measured perturbations. Others reason about controllability inside an inferred network. In our survey, none combined measured effects with this target-specific outside-the-cone certificate. That pairing—not the regulator list—is the novelty claim.” |
| 7 | 2:24–2:45 | Fig. 4, K562 transfer and doubles | “The same optimization code also runs on a different K562 CRISPRa screen and recovers a held-out CEBPA state at cosine point eight seven eight. But the doubles are the guardrail: measured combinations match the sum of singles at median cosine point seven one. Additivity helps; it is not truth.” |
| 8 | 2:45–3:00 | `app/index.html`, end title card | “The output is a decision for the next experiment: test a focused knockdown panel, add another modality, or seek a better dictionary. Not ‘this target works’—but ‘is this intervention class even pointing the right way?’” |

## Clean voiceover

> A GPS does not just suggest a route—it tells you when the road you need is not on the
> map. We built that second answer for cell engineering.
>
> Drug discovery pays for wrong turns late. This tool does not claim to predict clinical
> success. It asks a smaller question early: with the intervention effects we have actually
> measured, can we even point the cell in the direction we want?
>
> Start with a real genome-scale CRISPRi screen in primary human T cells: thirty-three
> thousand nine hundred eighty-three perturbation-by-condition profiles across more than
> ten thousand readout genes. Each perturbation leaves a measured fingerprint—some genes
> rise, others fall.
>
> We may combine those fingerprints, but we cannot reverse one and pretend we performed a
> negative knockdown. Under a first-order additivity model, all non-negative combinations
> form a cone. The target state is an arrow. The question is: how closely can this measured
> cone point along it?
>
> For the Th2-to-Th1 direction, the answer is: partly. The held-out cosine is point four four
> eight, higher than all sixty shuffled targets. Thirty-nine percent of target energy is
> captured by measured knockdown directions. Twenty-five percent is captured only by a
> sign-flipped activation proxy, and thirty-five percent by neither. The residual proves the
> full target arrow sits outside the measured knockdown cone.
>
> Its largest positive coordinates rank what the closest knockdown mixture still
> under-delivers. Those genes are not proven activation targets—they are a short list for
> CRISPRa or de-repression tests. Separately, the top greedy knockdown panels collapse to one
> hundred two unique candidates; forty-five lack a conventional drug handle and ten have a
> clinical-grade drug in the saved snapshot.
>
> We structured a survey of ninety-one prior methods. Some use measured perturbations.
> Others reason about controllability inside an inferred network. In our survey, none
> combined measured effects with this target-specific outside-the-cone certificate. That
> pairing—not the regulator list—is the novelty claim.
>
> The same optimization code also runs on a different K562 CRISPRa screen and recovers a
> held-out CEBPA state at cosine point eight seven eight. But the doubles are the guardrail:
> measured combinations match the sum of singles at median cosine point seven one.
> Additivity helps; it is not truth.
>
> The output is a decision for the next experiment: test a focused knockdown panel, add
> another modality, or seek a better dictionary. Not “this target works”—but “is this
> intervention class even pointing the right way?”

## Production checklist

1. Record the live HTML explorers at 1920×1080. Use the static figures only for the paper
   beats.
2. Keep three on-screen facts visible long enough to read: **0.448 held-out** · **above all
   60 shuffles** · **0.71 additivity median**.
3. Show “directional, screen-relative, additive model” once on screen. It prevents the word
   “reachable” from being heard as phenotypic rescue.
4. Put the dataset credit on the final card: Zhu et al. 2025; Marson/Pritchard labs; CZI
   Virtual Cells Platform.
5. Record voiceover first, then cut the screen capture to it.

## Number-consistency note

The flagship dense null is the only null number to narrate: held-out cosine **0.448**, above
all **60** shuffles, plus-one empirical **p = 1/61**; the corresponding **z ≈ 24** is a
descriptive standardized separation.

The atlas explorer also displays held-out cosine **0.446** and z **45** from an older
eight-shuffle screening run. The cosine difference is split/rounding variation. The z
difference does **not** scale lawfully with the number of shuffles; eight draws simply give
an unstable estimate of the null standard deviation. Keep the atlas value labelled
“eight-shuffle screening estimate,” and use the 60-shuffle result in narration.

## Optional 60-second cut

Use beats **1 → 3 → 4 → 8**: GPS hook, cone, flagship verdict, and next-experiment close.
