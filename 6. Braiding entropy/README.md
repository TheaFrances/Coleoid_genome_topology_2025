## Braiding calculation using Braidlab in MATLAB

This section describes how to generate and analyse chromosomal braiding maps based on macrosynteny alignments and sliding window genome comparisons. This enables the quantification of topological braiding patterns between aligned chromosomal regions, providing insights into the spatial organisation and conservation of 3D genome architecture.

### Prepare alignment files for braiding

To generate input `.map` files used for braiding, the script [`getHomChrWindowBraidMap.pl`](getHomChrWindowBraidMap.pl) was used. This script uses an input file of homologous chromosomes identified via any macrosynteny tool (example input file: [`eup-obi.mbh.psynt.tab`](eup-obi.mbh.psynt.tab), as well as a file of aligned genomes in tabular bed format (example input file: [`blastn.res.filt.all`](blastn.res.filt.all) made with any genome aligner (e.g. `megablast`)

Example command to generate the braiding input file:

```bash
perl getHomChrWindowBraidMap.pl eup-obi.mbh.psynt.tab blastn.res.filt.all 10000000 50 > esc-obi.10m.braid.map
```

This creates a map file with a 10 Mb window and 50 sliding window steps e.g. [esc-obi.10m.braid.map](../Test%20input%20files//esc-obi.10m.braid.map)

### Braiding analysis in MATLAB

This step requires MATLAB version 2025a or newer and [Braidlab toolbox](https://github.com/jeanluct/braidlab). Open MATLAB and navigate to the directory containing the `.map` file. Load and run [`braiding_commands.m`](braiding_commands.m) script in MATLAB, adjusting the file name to match the correct comparison and window size.



