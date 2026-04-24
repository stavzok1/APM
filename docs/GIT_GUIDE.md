# Git + GitHub workflow (APM)

This guide matches how you use Git locally, connect to GitHub with **SSH**, and move work between machines (laptop ↔ cluster).

---

## What `ssh-keygen` does (and what your command means)

```bash
ssh-keygen -t ed25519 -C "stav.zok@gmail.com"
```


| Part             | Meaning                                                                                                                                                                                                                                                |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `**ssh-keygen**` | Program that creates a **key pair** for SSH: one **private** key (secret, stays on your computer) and one **public** key (safe to upload).                                                                                                             |
| `**-t ed25519`** | **Type** of cryptography: **Ed25519** is a modern, short key format. GitHub recommends Ed25519 (or RSA with enough bits).                                                                                                                              |
| `**-C "…"`**     | **Comment** stored inside the public key. It does **not** authenticate you by itself; it is a **label** so you remember which machine or purpose the key is for. People often use their **email**; it appears at the end of the `.pub` line on GitHub. |


After running it you typically get:

- `~/.ssh/id_ed25519` — **private** key (never copy this to GitHub or email it).
- `~/.ssh/id_ed25519.pub` — **public** key (this is what you paste into **GitHub → Settings → SSH and GPG keys**).

**Why this connects you to GitHub:** GitHub trusts your computer when a connection is signed with the **private** key that matches a **public** key registered on your account. Passwords are not used for `git@github.com:…` URLs in the same way as HTTPS.

**Test:**

```bash
ssh -T git@github.com
```

You should see a message that mentions your GitHub username.

---

## Core ideas (read this once)


| Term                  | Meaning                                                                                                                                                                          |
| --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Repository (repo)** | Your project + its **history**. On disk it is your files **plus** the hidden `**.git/`** folder (all commits, branches, remotes).                                                |
| **Commit**            | A **snapshot** of the tracked files you **staged**, with a message and author. Commits form a chain (history).                                                                   |
| **Branch**            | A **movable pointer** to a commit (e.g. `main`). New commits move the branch forward.                                                                                            |
| `**HEAD`**            | “Where you are now” — usually the **latest commit on your current branch**.                                                                                                      |
| **Remote**            | A **nickname** for another copy of the repo (almost always **GitHub** for you). The usual nickname is `**origin`**.                                                              |
| `**origin`**          | **Not** a GitHub feature name — it is only a **convention**: the default name for “the main remote I push to / pull from.” You could call it `github`, but `origin` is standard. |
| **Push**              | Send your **new commits** on a branch **to** the remote (e.g. GitHub).                                                                                                           |
| **Pull**              | **Fetch** updates from the remote and **integrate** them into your current branch (download + merge, by default).                                                                |


**Important:** `git commit` only writes to `**.git/` on your machine`**. GitHub sees nothing until` **git push`**.

---

## One-time setup (identity + GitHub SSH)

### Who appears on commits (any email GitHub knows is fine)

```bash
git config --global user.name "Your Name"
git config --global user.email "your-email@example.com"
```

Git stores these **inside every commit** (author metadata). Use an address tied to your GitHub account if you want contributions to count on your profile (GitHub documents this under email settings).

### Create an SSH key (if you do not have one yet)

```bash
ssh-keygen -t ed25519 -C "your-email@example.com"
```

Press Enter to accept the default path (`~/.ssh/id_ed25519`) unless you want a separate key per project.

### Register the **public** key on GitHub

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the **single line** ending in your comment. GitHub: **Settings → SSH and GPG keys → New SSH key**.

### Create the empty repo on GitHub

On the website: **New repository** (e.g. `APM`). You can skip README if you already have local commits.

---

## Connect this project to GitHub (`origin`)

From your project root:

```bash
cd /home/stavz/masters/gdc/APM
```

### Add the remote (once per clone)

**SSH URL** (uses your SSH key; no HTTPS password):

```bash
git remote add origin git@github.com:YOUR_USERNAME/APM.git
```

- `**git remote add**` — “remember this URL under a short name.”
- `**origin**` — the short name (convention).
- `**git@github.com:USER/REPO.git**` — GitHub’s SSH address for that repo.

Check:

```bash
git remote -v
```

You should see `fetch` and `push` lines for `origin`.

### If you made a mistake

```bash
git remote remove origin
git remote add origin git@github.com:YOUR_USERNAME/APM.git
```

### First push (creates `main` on GitHub and sets tracking)

```bash
git push -u origin main
```


| Part                        | Meaning                                                                                                                                                    |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `**git push**`              | Upload commits.                                                                                                                                            |
| `**origin**`                | Push to the remote named `origin`.                                                                                                                         |
| `**main**`                  | Your **local** branch name to send.                                                                                                                        |
| `**-u`** (`--set-upstream`) | Remember: “local `**main`** tracks `**origin/main**`.” After this, on branch `main` you can run `**git push**` and `**git pull**` with no extra arguments. |


If GitHub’s default branch is not `main`, either rename locally (`git branch -m main`) or push with the branch name you use.

---

## Daily workflow: status → stage → commit

### See what changed

```bash
git status
```

- **Untracked** — new files Git is not watching yet (or ignored by `.gitignore`).
- **Modified** — tracked files changed.
- **Staged** — changes included in the **next commit**.

### Stage files (choose what goes in the next commit)

```bash
git add path/to/file.py           # one file
git add pipeline/                 # whole directory
git add -p                        # patch by patch (interactive)
```

### Save a snapshot

```bash
git commit -m "Short, clear description of what changed"
```

- `**git commit**` — creates a new commit from the **staged** changes only.
- `**-m "…"`** — commit message (required unless an editor opens).

### See history

```bash
git log --oneline --graph -15
```

---

## Branching: when and how

**Use a branch** when you want an experiment or feature isolated from `main` until it is ready.

### Create and switch

```bash
git switch -c feature-name    # create branch and switch to it (modern)
# older equivalent:
# git checkout -b feature-name
```

### List branches

```bash
git branch          # local
git branch -a       # local + remote-tracking
```

### Switch branch

```bash
git switch main
```

### Bring `main` up to date from GitHub, then merge your work

```bash
git switch main
git pull                          # get latest origin/main into main
git merge feature-name            # merge your branch into main
git push                          # upload updated main
```

(There is also `**git rebase**` for a linear history; learn it after you are comfortable with merge.)

### Delete a local branch (after merge)

```bash
git branch -d feature-name
```

---

## Pulling and pushing (after `-u` is set)

On branch `main` with upstream configured:

```bash
git pull          # download + integrate remote changes into current branch
git push          # upload your new commits on current branch
```

### If you never set upstream

```bash
git push -u origin main
```

### Get remote branches without merging yet

```bash
git fetch origin
```

Useful to see what exists on GitHub before you merge or switch.

---

## Typical scenarios (APM / cluster)


| Scenario                                 | What to do                                                                                          |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------- |
| **Finished a logical change**            | `git add …` → `git commit -m "…"` → `git push`                                                      |
| **Start work on cluster**                | `git clone git@github.com:USER/APM.git` → `cd APM` → copy/sync `**data/`** outside Git (e.g. rsync) |
| **Cluster: get latest code from laptop** | `git pull` on `main` (or your working branch)                                                       |
| **Laptop: get latest from cluster**      | If you only committed on cluster, `**git pull`** on laptop after cluster `**git push`**             |
| **Big files (HiChIP, VCF, etc.)**        | **Do not commit**; keep them under `data/` (ignored) or set `APM_WORKING_DIR` / sync with **rsync** |


---

## Safety checks before push

```bash
git status
git diff                          # unstaged changes
git diff --cached                 # staged changes
git log origin/main..HEAD         # commits you are about to push (after fetch)
```

---

## Quick reference table


| Command                   | Purpose                              |
| ------------------------- | ------------------------------------ |
| `git status`              | What changed / what is staged        |
| `git add <path>`          | Stage changes                        |
| `git commit -m "msg"`     | Create commit from staged changes    |
| `git log --oneline`       | History                              |
| `git remote -v`           | Show remotes (`origin` URLs)         |
| `git push -u origin main` | First push + set upstream            |
| `git push`                | Push current branch (needs upstream) |
| `git pull`                | Fetch + merge into current branch    |
| `git fetch`               | Download remote refs only            |
| `git switch -c name`      | New branch and switch                |
| `git switch main`         | Back to `main`                       |


---

## If something goes wrong

- **Pushed a secret by mistake:** rotate the secret; consider `git filter-repo` (advanced) — prevention (`.gitignore`, no tokens in repo) is easier.
- `**git push` rejected (non-fast-forward):** someone else updated GitHub; run `**git pull`** (resolve conflicts if any), then `**git push`** again.

For project-specific ignores and large paths, see `**.gitignore**` in the repo root.