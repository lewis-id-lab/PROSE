You are checking two independent extractions of the same results from a published paper against the paper itself. The paper is attached. The two extractions are labelled A and B. You are not told which one is the original, and you should not try to work it out.

For every item and every field listed, find the value in the paper and say which extraction matches it.

- `value_in_source`: the value as the paper prints it. Use null if the paper does not report it.
- `verdict`: "A" if only A matches the paper, "B" if only B matches, "both" if both match (for example the same value written differently, or the alleles named the other way round with everything else oriented to match), "neither" if neither matches, and "not_found" if you cannot locate the value in the paper.
- `quote`: the exact text or table row you read the value from, copied character for character. Required unless the verdict is not_found.
- `page_or_table`: where you read it.

Judge only against the paper. Do not prefer an extraction because it looks more complete or more plausible, and do not recompute anything from counts. If the two extractions describe different results (for example different subgroups or different genetic models), say so in `item_note` and judge each field against the result that the item's `paper_result` describes.
