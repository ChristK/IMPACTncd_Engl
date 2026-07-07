# Running IMPACTncd England on your own computer — a beginner's guide

This guide takes you from a **fresh Windows, macOS, or Linux machine** all the
way to **running your own simulations** with the IMPACTncd England model. No
programming background is assumed. Everything the model needs — the R
environment, ~13 GB of input data, and the model code — is bundled inside a
single **Docker image**, so you do not install R or compile anything by hand.

> **In one sentence:** install Docker → make a folder with a launcher script and
> a config file → run one command. The launcher downloads the model image and
> drops you (or your script) inside a ready-to-run container.

If you are on the project's shared Liverpool server instead of your own machine,
use the server handout (`TEAM_QUICKSTART`) rather than this guide.

---

## Contents

1. [Before you start — what you need](#1-before-you-start--what-you-need)
2. [Install Docker](#2-install-docker)
3. [Give Docker enough memory and disk](#3-give-docker-enough-memory-and-disk)
4. [Check Docker works](#4-check-docker-works)
5. [Make your project folder and get the launcher](#5-make-your-project-folder-and-get-the-launcher)
6. [Edit your configuration file](#6-edit-your-configuration-file)
7. [Write a run script](#7-write-a-run-script)
8. [Run the model](#8-run-the-model)
9. [Find your results](#9-find-your-results)
10. [Run a policy scenario](#10-run-a-policy-scenario)
11. [Read the built-in tutorials (vignettes)](#11-read-the-built-in-tutorials-vignettes)
12. [Troubleshooting](#12-troubleshooting)
13. [Updating and cleaning up](#13-updating-and-cleaning-up)
14. [How it works (optional background)](#14-how-it-works-optional-background)

---

## 1. Before you start — what you need

How much hardware you need depends entirely on **what you are doing**:

- **Just trying it out** — a handful of iterations, to confirm the model runs on
  your machine. A typical laptop is fine.
- **Real analysis** — a meaningful, stable result needs **hundreds of
  Monte-Carlo iterations**. That is a far heavier job, and hardware matters a
  lot (see the warning below).

| Requirement | Just trying it out | Real analysis |
|---|---|---|
| **Operating system** | Windows 10/11 (64-bit), macOS (Intel or Apple Silicon), or Linux | same |
| **Disk space** | **~40 GB free** — the image is ~33 GB unpacked, plus room for your results and the synthetic-population cache | same, plus more for larger result sets |
| **CPU cores** | 1 core works | **as many as you can** — the run is split across cores, so more cores finish proportionally sooner |
| **Memory (RAM)** | ~10–16 GB | **plenty** — the simulation uses roughly **10 GB per core**, so more cores means more RAM (e.g. 8 cores ≈ 80 GB) |
| **Time** | first run: several minutes to a couple of hours (it builds the synthetic population once); a 2-iteration test after that: minutes | **hours to days** even on a well-specced machine — see the warning below |
| **Internet** | once, to download Docker and the image | same |

> ### ⚠️ Realistic runtimes — read this before planning real work
>
> A publishable result typically needs **hundreds of Monte-Carlo iterations**.
> On a modest setup — **16 GB of RAM and a single CPU core** — that many
> iterations can take **days or even weeks** to finish. The simulation
> parallelises well, so the practical answer is hardware: **more CPU cores and
> plenty of RAM give far more reasonable runtimes and a much smoother
> experience.** Because each core needs ~10 GB of RAM, scaling up the number of
> cores means scaling up memory too. If all you have is a laptop, use it to
> learn the workflow with a couple of iterations, and run the large jobs on a
> workstation, a shared server, or a cloud VM.

You do **not** need a GitHub account, a Zenodo account, or any access token —
everything used here is public.

---

## 2. Install Docker

Docker is the program that runs the packaged model. Install it for your OS.

### Windows

1. Download and install **Docker Desktop for Windows**:
   <https://docs.docker.com/desktop/install/windows-install/>
2. When asked, keep the **WSL 2 based engine** option ticked (recommended).
   If Windows prompts you to install/update WSL, accept it and reboot.
3. Start **Docker Desktop** from the Start menu and wait until the whale icon in
   the system tray stops animating (it says *"Docker Desktop is running"*).

You will type commands into **PowerShell** (search "PowerShell" in the Start
menu). PowerShell 5.1 comes with Windows 10/11 — nothing to install.

### macOS

1. Download and install **Docker Desktop for Mac** (pick the Apple Silicon or
   Intel build to match your Mac):
   <https://docs.docker.com/desktop/install/mac-install/>
2. Start **Docker Desktop** from Applications and wait until the whale icon in
   the menu bar shows *"Docker Desktop is running"*.
3. Install a small helper the launcher script needs (`realpath`), using
   [Homebrew](https://brew.sh):
   ```bash
   brew install coreutils
   ```
   You will type commands into the **Terminal** app (Applications → Utilities → Terminal).

### Linux

Install **Docker Engine** using your distribution's official instructions:

- [Ubuntu](https://docs.docker.com/engine/install/ubuntu/) ·
  [Debian](https://docs.docker.com/engine/install/debian/) ·
  [Fedora](https://docs.docker.com/engine/install/fedora/) ·
  [CentOS / Rocky Linux](https://docs.docker.com/engine/install/centos/)

Then let your user run Docker without `sudo`:

```bash
sudo usermod -aG docker $USER
```

**Log out and back in** (or reboot) for that group change to take effect.

---

## 3. Give Docker enough memory and disk

The model is memory-hungry — it uses roughly **10 GB of RAM per CPU core**. If
Docker isn't allowed enough memory, the simulation's workers are killed mid-run
and you get a confusing downstream error (see the entry in
[Troubleshooting](#12-troubleshooting)). **How you raise the memory depends on
your OS — and, on Windows, on which Docker backend you use.**

**Disk (all platforms):** make sure there is **~40 GB** free for the image plus
your results.

### macOS

Docker Desktop → **Settings (gear icon) → Resources → Memory**: set it to at
least **12 GB** (more if you have it), then **Apply & restart**.

### Windows

Windows Docker Desktop has **two backends** with *different* memory controls.
Check which one you have at **Settings → General →** the *"Use the WSL 2 based
engine"* checkbox:

- **WSL 2 based engine — ticked (this is the default).** There is **no memory
  slider** in Docker Desktop; memory is managed by WSL 2, which by default caps
  itself at only **~50% of your RAM or 8 GB, whichever is smaller**. On a 16 GB
  laptop that's an 8 GB ceiling — *below* the ~10 GB one core needs, which is the
  usual reason runs die. Raise it by creating the file
  `C:\Users\<your-username>\.wslconfig` (Notepad is fine) containing:
  ```ini
  [wsl2]
  memory=12GB
  ```
  Then, in PowerShell, apply it and restart Docker Desktop:
  ```powershell
  wsl --shutdown
  ```
  (Use more than 12 GB if your machine has it, e.g. `memory=24GB`. Leave a few GB
  for Windows itself.)

- **Hyper-V engine — unticked (legacy).** Here the slider *does* work: Docker
  Desktop → **Settings → Resources → Memory** → at least **12 GB** →
  **Apply & restart**.

### Linux

Docker uses your machine's memory and disk directly — there is no separate
setting. Just make sure you have enough free RAM and disk (see
[section 1](#1-before-you-start--what-you-need)).

---

## 4. Check Docker works

Open your terminal (PowerShell on Windows; Terminal on macOS/Linux) and run:

```bash
docker info
```

If Docker is running you get a page of information and **no error**. If instead
you see *"Cannot connect to the Docker daemon"*:

- **Windows / macOS:** make sure Docker Desktop is actually started (whale icon).
- **Linux:** start it with `sudo systemctl start docker`, and confirm you did
  the `usermod` step and logged out/in.

---

## 5. Make your project folder and get the launcher

Everything for one project lives in a single folder that holds **three files**:
the launcher script, a configuration file, and your run script. Create the
folder, move into it, and download the launcher plus a starter configuration.

> **Tip:** paste these **one line at a time** — paste a line, press Enter, wait,
> then the next. Pasting a whole block at once can jam commands together in some
> terminals.

### Linux / macOS (Terminal)

```bash
mkdir -p ~/impactncd_run
cd ~/impactncd_run
curl -O https://raw.githubusercontent.com/ChristK/IMPACTncd_Engl/main/docker_setup/setup_user_docker_env.sh
curl -O https://raw.githubusercontent.com/ChristK/IMPACTncd_Engl/main/inputs/sim_design.yaml
```

### Windows (PowerShell)

```powershell
mkdir $HOME\impactncd_run
cd $HOME\impactncd_run
Invoke-WebRequest https://raw.githubusercontent.com/ChristK/IMPACTncd_Engl/main/docker_setup/setup_user_docker_env.ps1 -OutFile setup_user_docker_env.ps1
Invoke-WebRequest https://raw.githubusercontent.com/ChristK/IMPACTncd_Engl/main/inputs/sim_design.yaml -OutFile sim_design.yaml
```

You now have a folder containing `setup_user_docker_env.sh` (or `.ps1` on
Windows) and `sim_design.yaml`. There is **no need to make the script
executable** — you will run it with `bash` (or PowerShell), which does not care.

> **Keep only one `sim_design*.yaml` in this folder.** The launcher finds your
> configuration automatically by that name; if there are several, it will ask
> you to pick one.

---

## 6. Edit your configuration file

`sim_design.yaml` controls the simulation. Open it in any text editor (Notepad,
TextEdit, VS Code, `nano` …). The file has many options; for a first run you
only need to change **four lines**. Find them and set them like this:

```yaml
output_dir: ./outputs
synthpop_dir: ./synthpop
clusternumber: 1
clusternumber_export: 1
```

What these mean:

- **`output_dir`** — where results are written. `./outputs` means a folder named
  `outputs` **inside your project folder** — you do not create it, the launcher
  does. Your results end up right next to your config.
- **`synthpop_dir`** — where the synthetic population is cached (built on the
  first run, reused afterwards). `./synthpop` keeps it in your project folder too.
- **`clusternumber`** — how many CPU cores the simulation uses. Each core needs
  **~10 GB of RAM**, so `1` is the safe choice on a personal machine. Raise it
  only if you have plenty of RAM (e.g. `2` needs ~20 GB).
- **`clusternumber_export`** — the same idea for the results-export step. Keep it
  at `1` as well; the default is high and would try to use a lot of memory.

> The downloaded file ships with the model author's own paths and large core
> counts — that is why changing these four lines matters. Leave every other
> setting at its default for now.

Save the file.

---

## 7. Write a run script

Create a new file called **`simulate.R`** in the same folder, with this content
(a minimal baseline run):

```r
# simulate.R — a minimal baseline run
source("global.R")                          # load the model (this file is inside the image)

IMPACTncd <- Simulation$new("sim_design.yaml")   # read the config next to this script

IMPACTncd$run(1:2, multicore = TRUE, "sc0")      # 2 iterations of the baseline "sc0"
IMPACTncd$export_summaries(multicore = TRUE)     # write summary tables

cat("Done! Results are in your ./outputs folder.\n")
```

Your project folder now looks like this:

```text
impactncd_run/
  setup_user_docker_env.sh     (or .ps1 on Windows)   the launcher
  sim_design.yaml                                      your configuration
  simulate.R                                           your run script
```

> **Why `1:2`?** That runs just two Monte-Carlo iterations — enough to confirm
> everything works quickly. **Real analysis needs hundreds** (e.g. `1:200`) for
> stable estimates — which is exactly why hardware matters so much; see the
> [runtime warning in section 1](#1-before-you-start--what-you-need).

---

## 8. Run the model

Run the launcher **from inside your project folder**. It auto-detects your
`sim_design.yaml`, downloads the model image the first time (be patient), and
then either drops you into an interactive session or runs your script for you.

### Linux / macOS

```bash
cd ~/impactncd_run
bash setup_user_docker_env.sh -Run
```

### Windows (PowerShell)

```powershell
cd $HOME\impactncd_run
# allow scripts to run in THIS PowerShell window only (needed once per window):
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\setup_user_docker_env.ps1 -Run
```

**What `-Run` does:** it runs the single `simulate.R` in your folder from start
to finish and then exits — hands-off. When it finishes you are back at your
normal prompt and your results are on disk.

### Prefer to work interactively?

Leave off `-Run` (`bash setup_user_docker_env.sh`, or
`.\setup_user_docker_env.ps1`). Instead of running your script, the launcher
drops you into a shell **inside** the container. Your prompt changes to something
like `youruser@a1b2c3d4:/IMPACTncd_England$`. Your files (`sim_design.yaml`,
`simulate.R`) are right there next to the model — no subfolder to remember — so
you can run:

```r
R                                       # start R
source("global.R")
IMPACTncd <- Simulation$new("sim_design.yaml")
IMPACTncd$run(1:2, multicore = TRUE, "sc0")$export_summaries(multicore = TRUE)
```

Type `q(save = "no")` to leave R, then `exit` to leave the container.

> **First-run messages are normal.** You will see lines about pulling the image,
> `Auto-detected sim-design YAML …`, `Mounting output_dir …`, a few
> `Entrypoint debug` lines, and a note that inside the container the folders are
> called `/outputs` and `/synthpop`. All expected.

> **Windows/macOS performance tip:** if runs feel very slow, add `-UseVolumes`
> (PowerShell) or `--UseVolumes` (bash). This uses faster Docker-managed storage
> and copies your results back to `./outputs` when the container exits.

---

## 9. Find your results

Results are written to the **`outputs`** folder inside your project folder — on
your real machine, not lost inside the container:

```text
impactncd_run/outputs/
  summaries/        tables such as ...  (*.csv.gz)
  ...
```

Open the `.csv.gz` files with R, Python, or a spreadsheet that reads gzipped
CSV. The `synthpop` folder holds the cached synthetic population — leave it in
place so future runs are fast; delete it only if you want to rebuild from scratch.

> **To redo a run cleanly**, empty the `outputs` folder first (or add
> `IMPACTncd$del_logs()$del_outputs()` before the `run(...)` line in
> `simulate.R`) — otherwise leftovers from the previous attempt mix into the new
> results.

---

## 10. Run a policy scenario

The point of the model is to compare a **baseline** with a **what-if policy**.
Replace `simulate.R` with the version below. It runs the baseline (`sc0`), then
applies a policy (10 % lower BMI from 2023 onward) and runs that (`sc1`):

```r
# simulate.R — baseline vs. a BMI policy
source("global.R")
IMPACTncd <- Simulation$new("sim_design.yaml")

# 1) Baseline
IMPACTncd$run(1:2, multicore = TRUE, "sc0")

# 2) Policy: 10% lower BMI from year 2023 onward
IMPACTncd$update_primary_prevention_scn(
  function(sp) sp$pop[year >= 23L, bmi_curr_xps := bmi_curr_xps * 0.9]
)
IMPACTncd$run(1:2, multicore = TRUE, "sc1")

IMPACTncd$export_summaries(multicore = TRUE)
cat("Done! Compare sc0 vs sc1 in ./outputs.\n")
```

Run it exactly as in [section 8](#8-run-the-model). Other risk factors you can
nudge the same way include blood pressure, cholesterol, smoking, alcohol, fruit
and vegetable intake, and physical activity. For the full menu, read the
built-in tutorials — see the [next section](#11-read-the-built-in-tutorials-vignettes).

---

## 11. Read the built-in tutorials (vignettes)

The model ships several tutorials ("vignettes"). List them from inside R with:

```r
vignette(package = "IMPACTncdEngland")
```

Note that **`vignette("how_to_test_run")` will fail** inside the container with
`sh: 1: xdg-open: not found`. That is expected — the container has no web browser
to open the page in. Nothing is broken; you view the tutorials on your own
computer instead, like this:

**Copy the tutorials into your results folder, then open them in your browser.**
In your R session inside the container:

```r
doc <- system.file("doc", package = "IMPACTncdEngland")
dir.create("/outputs/vignettes", showWarnings = FALSE)
file.copy(list.files(doc, full.names = TRUE), "/outputs/vignettes", overwrite = TRUE)
```

Then, **on your own computer**, open your project folder → `outputs/vignettes/`
and double-click **`index.html`** — it links to every tutorial, nicely
formatted. (`/outputs` inside the container is the very same folder as
`./outputs` on your machine, so the copied files are already on your disk.)

The tutorials available:

| Vignette | What it covers |
|---|---|
| `how_to_test_run` | **Start here** — running a test simulation. |
| `how_to_run_scenarios` | Defining and running a baseline plus policy scenarios. |
| `understanding_model_outputs` | Interpreting the output files. |
| `custom-scenario-columns` | Exporting custom columns created in scenarios. |
| `inputs_manifest_system` | The inputs manifest / data-asset tracking system. |
| `zenodo_data_management` | Downloading and managing the model's input data. |

> **Prefer to stay in the terminal?** The tutorial sources are plain text you can
> read in place — for example
> `less /usr/local/lib/R/site-library/IMPACTncdEngland/doc/how_to_test_run.Rmd`
> (press `q` to quit). The matching `.R` file beside it holds just the runnable
> code from that tutorial.

---

## 12. Troubleshooting

**"Cannot connect to the Docker daemon."**
Docker is not running or you lack permission. Start Docker Desktop
(Windows/macOS); on Linux run `sudo systemctl start docker` and make sure you
did the `usermod -aG docker $USER` step and logged out/in.

**Windows: "running scripts is disabled on this system."**
PowerShell's execution policy blocked the `.ps1`. Run this in the same window,
then re-run the launcher:
```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

**macOS: `realpath: command not found` (or the script misbehaves early).**
Install coreutils: `brew install coreutils`.

**"Failed to pull Docker image."**
Check your internet connection and that Docker is running. The image is
`chriskypri/impactncdengl:main` — you can confirm it exists at
<https://hub.docker.com/r/chriskypri/impactncdengl/tags>.

**The run fails with `... did not deliver results` and/or `Invalid Input Error:
Provided table/dataframe must have at least one column` (on `SELECT DISTINCT mc
FROM lc_table`) — or your machine simply freezes.**
The simulation's parallel workers were **killed mid-run, almost always out of
memory**, so no results were written and the export step then chokes on the
empty dataset. Fix it in this order:
1. **Give Docker more memory** — see
   [section 3](#3-give-docker-enough-memory-and-disk). On Windows/WSL 2 the
   default cap is often just 8 GB, below the ~10 GB one core needs; raise it via
   `.wslconfig`.
2. **Keep `clusternumber` and `clusternumber_export` at `1`** in your YAML — each
   extra core needs another ~10 GB.
3. **Pinpoint the cause** by running single-threaded: set `multicore = FALSE` in
   both your `run(...)` and `export_summaries(...)` calls. That removes the
   parallel layer, so a real error (if there is one) is printed directly instead
   of the vague "did not deliver results". If it then succeeds, the problem was
   memory/cores — go back to steps 1–2.

**The first run takes a very long time.**
Expected. The first run downloads the image and builds the synthetic population
once. Subsequent runs reuse the cached `synthpop` folder and are much faster.

**`vignette(...)` gives `xdg-open: not found`.**
Expected — the container has no browser to open the page in. See
[section 11](#11-read-the-built-in-tutorials-vignettes) for how to read the
tutorials on your own computer instead.

**"multiple YAML files … cannot choose one automatically."**
You have more than one `sim_design*.yaml` in the folder. Remove the extras, or
add `-SimDesignYaml sim_design.yaml` to the launch command to pick one.

**Windows: `WARNING: wsl.exe not found; falling back to ... '/c/...'`.**
The launcher couldn't run `wsl.exe` to translate your Windows path, so it used
the legacy `/c/...` form. If you use the **WSL 2 backend** (the default), that
form may not map to your real folder, so your `outputs/` can end up empty on the
host *even when the run succeeds*. Make sure WSL is installed and `wsl.exe` is on
your PATH (it lives at `C:\Windows\System32\wsl.exe` — a normal PowerShell window
finds it). If you deliberately use the Hyper-V backend, the `/c/...` form is
correct and this warning is harmless.

---

## 13. Updating and cleaning up

**Get the latest model.** The image is updated over time. To fetch the newest
version, pull it (the launcher also does this automatically each run):
```bash
docker pull chriskypri/impactncdengl:main
```
It is also good practice to re-download the launcher script occasionally so you
have the latest features (repeat the `curl`/`Invoke-WebRequest` step from
[section 5](#5-make-your-project-folder-and-get-the-launcher)).

**Free up disk space.** The image is large (~33 GB). Remove it when you are done:
```bash
docker rmi chriskypri/impactncdengl:main
```
Remove **all** unused Docker data (⚠️ affects every image/container on the
machine):
```bash
docker system prune -a
```
Your results and config are ordinary files in your project folder — deleting the
image does not touch them.

---

## 14. How it works (optional background)

- **Image vs. container.** The *image* is a frozen, read-only copy of the whole
  model (software + data). Each time you run the launcher it starts a fresh,
  disposable *container* from that image. The container vanishes when it exits —
  only your mounted folders (`outputs`, `synthpop`) survive, because they live on
  your real disk.
- **Why your files "appear" next to the model.** The launcher mounts your project
  folder into the container and links your files into the model's own folder, so
  you can refer to `sim_design.yaml` and `simulate.R` with no path prefix. If one
  of your files happened to share a name with a model file, yours is skipped (the
  model's copy wins) and a note is printed.
- **Runs as *you*.** Files the container creates (your results) are owned by your
  user account, not by `root`, so you can read and delete them normally.
- **Where the data comes from.** The ~13 GB of inputs and pre-computed tables are
  published on Zenodo (concept DOI `10.5281/zenodo.20812409`, CC-BY-SA-4.0) and
  are already baked into the image — you never download them by hand.

**Prefer a native R install instead of Docker?** That path (install R, build the
package, download the data from Zenodo) is documented in the main project
[`README.md`](../README.md) under *Option B — Native R installation*.

For the full reference on the launcher's options, the three-layer image
architecture, and building images yourself, see
[`docker_setup/README.md`](README.md).

---

*Questions or problems? Open an issue at
<https://github.com/ChristK/IMPACTncd_Engl/issues>.*
