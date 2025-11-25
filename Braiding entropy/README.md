## Braiding calculation using Braidlab in MATLAB

This section describes how to generate and analyse chromosomal braiding maps based on macrosynteny alignments and sliding window genome comparisons. This enables the quantification of topological braiding patterns between aligned chromosomal regions, providing insights into the spatial organisation and conservation of 3D genome architecture.

### Requirements

* MATLAB version 2025a or newer
* [Braidlab toolbox](https://github.com/jeanluct/braidlab)

### Step 1: Prepare alignment files for braiding

To generate input `.map` files used for braiding:

* Identify homologous chromosomes between two species using a macrosynteny tool (e.g. output file: `eup-obi.mbh.psynt.tab`)
* Align genomes using a whole-genome aligner such as `megablast`, and store aligned regions in tabular BED-like format (e.g. `blastn.res.filt.all`)

Generate the braiding input file using:

```bash
perl getHomChrWindowBraidMap.pl eup-obi.mbh.psynt.tab blastn.res.filt.all 10000000 50 > esc-obi.10m.braid.map
```

This creates a map file with a 10 Mb window and 50 sliding window steps.

### Step 2: Braiding analysis in MATLAB

* Open MATLAB
* Navigate to the directory containing the `.map` file
* Load and run [`braiding_commands.mat`](braiding_commands.mat) script in MATLAB, adjusting the file name to match the correct comparison and window size.

Example usage:

```matlab
load('esc-obi.10m.braid.map');
% Run custom script or Braidlab functions to analyze braiding patterns
```

---


