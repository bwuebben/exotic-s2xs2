# Getting, building, and installing GAP (macOS)

Recipe used on this machine 2026-07-09 (GAP 4.16.0, Apple Silicon, macOS 26). Kept here
because the answer to "how do I install GAP on a Mac" is genuinely non-obvious.

## The situation (verified 2026-07-09)

- **Homepage:** https://www.gap-system.org — releases at https://github.com/gap-system/gap/releases
- **Homebrew: no formula or cask exists** in homebrew-core, under any name. `brew install gap` will never work.
- **MacPorts: no port** (ports.macports.org API returns zero results).
- The **official Homebrew tap** `gap-system/gap` exists (`brew install gap-system/gap/gap`)
  but ships **no bottles** — it source-compiles the *full* ~566 MB distribution plus heavy
  deps (Singular, PARI, fplll, zeromq, …) and runs `BuildPackages.sh` over all ~150 bundled
  packages (20–40 min, occasionally flaky). Use only if you want brew to manage upgrades.
- **GAP publishes no macOS binary.** The only prebuilt Mac binary anywhere is **Gap.app**
  (https://cocoagap.sourceforge.io/), a GUI frontend — last release Nov 2022, bundling the
  now-stale GAP 4.12.1. **conda-forge** has native arm64 binaries (`gap-defaults`) if you
  already run conda; there is also a Docker image (`gapsystem/gap-docker`, amd64).
- On **Linux** none of this matters: `apt install gap` / `dnf install gap` just works.

So on a clean Mac the right move is the **minimal source build** below: core tarball +
the four required packages. ~130 MB of downloads, ~2 minutes of compile on an M-series
CPU, and it is the *latest* GAP rather than the tap's N−1 or Gap.app's 2022 vintage.

## Minimal source build

Prerequisites: Xcode command-line tools, plus GNU gmp and readline from brew
(GAP cannot use macOS's libedit):

```bash
xcode-select --install          # if not already present
brew install gmp readline
```

Fetch the latest release tag, then download the **core** tarball and the **required
packages** tarball (gapdoc, primgrp, smallgrp, transgrp — GAP will not start without
them) and verify checksums:

```bash
V=$(curl -s https://api.github.com/repos/gap-system/gap/releases/latest | sed -n 's/.*"tag_name": "v\(.*\)".*/\1/p')
mkdir -p ~/opt/src && cd ~/opt/src
for f in gap-$V-core.tar.gz packages-required-v$V.tar.gz; do
  curl -sLO https://github.com/gap-system/gap/releases/download/v$V/$f \
       -sLO https://github.com/gap-system/gap/releases/download/v$V/$f.sha256
  echo "$(awk '{print $1}' $f.sha256)  $f" | shasum -a 256 -c
done
```

Extract — **gotcha:** the packages tarball has the package dirs at its *top level*; it
must be extracted *inside* `pkg/`, not at the GAP root:

```bash
cd ~/opt && tar xzf src/gap-$V-core.tar.gz
cd gap-$V && mkdir -p pkg && tar xzf ../src/packages-required-v$V.tar.gz -C pkg
```

Build (point configure at brew's gmp/readline) and put `gap` on PATH:

```bash
./configure --with-gmp=/opt/homebrew/opt/gmp --with-readline=/opt/homebrew/opt/readline
make -j"$(sysctl -n hw.ncpu)"        # ~2 min on Apple Silicon
ln -sf ~/opt/gap-$V/gap /opt/homebrew/bin/gap
```

No `make install` — GAP runs happily out of its build directory; the symlink is enough.

Sanity check:

```bash
echo 'Print(GAPInfo.Version, " ", StructureDescription(SmallGroup(8,3)), "\n");' | gap -q -A
# expect:  4.16.0 D8
```

## Notes

- The core library already contains everything this project's scripts use (finitely
  presented groups, Todd–Coxeter coset enumeration, `AbelianInvariants`, `GQuotients`,
  `PSL`/`SL`/`AlternatingGroup`, `RandomSource`). No optional packages needed.
- Run scripts non-interactively with `gap -q -A script.g` (`-q` quiet, `-A` skip
  autoloading optional packages); scripts here end with `QUIT_GAP(0);`.
- Reproducibility: `RandomSource(IsMersenneTwister, <seed>)` streams are stable across
  GAP versions — the memo's seeded 120-trial experiment reproduced exactly on 4.16.0.
