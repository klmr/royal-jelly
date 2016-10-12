## Royal jelly sequencing project

Honeybees feed larvae a royal jelly diet to trigger their development into queen
bee adults. The mechanism behind this nutrient-based differential development
[has been hypothesised][1] to be caused by Major Royal Jelly Protein 1 (MRJP1,
aka. royalactin). This is currently the focus of a scientific controversy (see
[2][], [3][]).

Here we investigate the presence of potentially regulatory RNAs (small RNAs and
long RNAs) in royal jelly, which might likewise have an effect on differential
development.

## Prerequisites

* GNU Make 3.81
* R 3.3.1

## Usage

The analysis pipeline is controlled via a Makefile. A list of commands can be
displayed via

```bash
make
```

All generated output is stored under `data`. In particular, `data/qc` contains
QC reports.

[1]: http://europepmc.org/abstract/MED/21516106 "Kamakura M., Royalactin induces queen differentiation in honeybees (2011)"
[2]: http://europepmc.org/abstract/MED/27652566 "Buttstedt A. & al., Royalactin is not a royal making of a queen (2016)"
[3]: http://europepmc.org/abstract/MED/27652567 "Kamakura M., Kamakura replies (2016)"
