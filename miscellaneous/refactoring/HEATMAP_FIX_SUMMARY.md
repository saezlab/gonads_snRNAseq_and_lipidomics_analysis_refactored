# Heatmap PDF Corruption Fix

## Problem
The workflow was generating corrupted PDF files when creating pathway enrichment heatmaps at line 488 of `workflows/02_fcg_pseudobulk_analysis.R`.

## Root Cause
**Namespace conflict between `pheatmap` and `ComplexHeatmap` packages.**

In `scripts/visualization/heatmaps.R`, both packages were loaded:
```r
library(pheatmap)        # Line 13
library(ComplexHeatmap)  # Line 14
```

When both packages are loaded, `ComplexHeatmap::pheatmap()` **masks** the original `pheatmap::pheatmap()` function. The ComplexHeatmap version:
1. Returns a different object type (Heatmap instead of gtable)
2. Doesn't properly handle the `filename` parameter
3. Requires explicit drawing with `draw()` when saving to PDF
4. Was causing PDF files to be created but not properly written, resulting in corruption

## Solution
**Use explicit namespace qualification: `pheatmap::pheatmap()`**

This ensures we always call the original pheatmap function, regardless of which packages are loaded in the namespace.

### Files Modified

#### 1. `scripts/visualization/heatmaps.R`

**Changes:**
- **Line 153**: `create_expression_heatmap()` now uses `pheatmap::pheatmap()`
- **Line 261**: `create_enrichment_heatmap()` now uses `pheatmap::pheatmap()`
  - Removed complex branching logic for ComplexHeatmap vs pheatmap
  - Simplified to always use original pheatmap with explicit namespace
  - Added file size validation after creation
  - Improved error handling with corrupted file cleanup
- **Line 356**: `create_correlation_heatmap()` now uses `pheatmap::pheatmap()`

**Key improvement in `create_enrichment_heatmap()`:**
```r
# OLD CODE (BROKEN):
if (using_complexheatmap) {
  pdf(filename, ...)
  p <- pheatmap(...)  # Could be ComplexHeatmap version!
  dev.off()
}

# NEW CODE (FIXED):
p <- pheatmap::pheatmap(  # Always uses original pheatmap
  enrichment_matrix,
  filename = filename,  # Handles PDF creation internally
  ...
)
# Validates file was created successfully
```

#### 2. `workflows/02_fcg_pseudobulk_analysis.R`

**Additional improvements made:**
- Added `tryCatch()` error handling around heatmap creation (lines 468, 543)
- Added matrix validation before plotting (lines 476, 551)
- Improved error messages to identify which heatmap fails

## Testing

### Test Script
Run `test_heatmap_fix_final.R` to verify the fix:

```r
source("test_heatmap_fix_final.R")
```

This will:
1. Load existing pathway results
2. Create a test heatmap PDF
3. Validate the PDF file is not corrupted (size > 5KB)
4. Report success/failure

### Expected Output
```
✓ Heatmap function completed without error
✓ File size looks good (> 5KB)
✓ PDF is likely NOT corrupted

File: results/pseudobulk/pathways/TEST_pathway_enrichment_heatmap.pdf
```

### Verification
Open the generated PDF file:
- **If it opens and displays the heatmap correctly**: Fix is successful! ✓
- **If it fails to open or is blank**: There may be additional issues

## Why This Happened

The issue was subtle because:
1. The function appeared to run without errors
2. A PDF file was created (so file.exists() returned TRUE)
3. But the PDF content was incomplete/corrupted
4. The "if (totalCols == 0)" error was coming from ComplexHeatmap's internal checks

## Prevention

**Best Practice**: When using multiple packages with overlapping function names, always use explicit namespace qualification:

```r
# Good - explicit
pheatmap::pheatmap(data, ...)
ComplexHeatmap::Heatmap(data, ...)

# Bad - ambiguous
pheatmap(data, ...)  # Which package's function?
```

## Additional Notes

- The workflow now properly handles both pathway and TF activity heatmaps
- All corrupted files are automatically deleted if creation fails
- File size validation helps catch corruption early
- The fix maintains backward compatibility - same function signatures

## Files to Test

After running workflow 02, verify these files are created and not corrupted:

1. `results/pseudobulk/pathways/pathway_enrichment_heatmap.pdf` (~8-12 KB)
2. `results/pseudobulk/pathways/tf_activity_heatmap.pdf` (~3-5 KB)

Both should open correctly in a PDF viewer and display colorful heatmaps.
