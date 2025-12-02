# 🤝 Contributing to Cecret

We welcome contributions to the Cecret Nextflow workflow from the community! By participating in this project, you agree to abide by our Code of Conduct.

## ❓ How to Contribute

There are several ways you can help improve Cecret:

1.  **Report a Bug:** Find an issue or error in the pipeline.
2.  **Suggest an Enhancement:** Propose a new feature, a new organism workflow, or an updated tool version.
3.  **Improve Documentation:** Fix typos, clarify usage, or add new tutorials.
4.  **Submit Code Changes (Pull Request):** Write and submit code to fix a bug or implement a feature.

---

## 🐛 Reporting Bugs

Before reporting a bug, please search [existing issues](https://github.com/UPHL-BioNGS/Cecret/issues) to see if the problem has already been addressed.

When submitting a new bug report via a GitHub Issue, please include the following critical information:

* **Cecret Version:** The specific Git tag or commit hash you are using (e.g., `v1.5.0` or `main` branch).
* **Execution Profile:** The Nextflow profile used (e.g., `-profile docker`, `-profile singularity`).
* **Operating System/Environment:** (e.g., Ubuntu 22.04, HPC Slurm cluster).
* **The Exact Error:** Copy and paste the complete traceback or error message.
* **Nextflow Log:** Attach or paste the content of the `.nextflow.log` file.
* **Reproduction Steps:** Provide the minimal `nextflow run` command and any necessary input files (or links to public accessions) that can reliably recreate the bug.

---

## ✨ Suggesting Enhancements

If you have an idea for a new feature or an improvement to an existing process:

* Use the [GitHub Issues page](https://github.com/UPHL-BioNGS/Cecret/issues) and clearly label the issue as a **"Feature Request."**
* Describe the goal of the feature and why it would be beneficial to the Cecret community.
* Outline how the feature should be implemented, particularly if it involves adding new parameters, tools, or configuration files.

---

## 💻 Code Contribution Guidelines

We use the standard GitHub **fork and pull request** workflow.

### Setup and Workflow

1.  **Fork** the repository to your personal GitHub account.
2.  **Clone** your fork locally: `git clone git@github.com:YOUR-USERNAME/Cecret.git`
3.  **Create a dedicated branch** for your contribution: `git checkout -b fix/issue-123-bug` or `feature/add-new-tool`.
4.  **Install dependencies** and run a minimal test of the pipeline locally before coding to ensure your environment is set up correctly.

### Coding Standards

* **Nextflow/Groovy:** Follow standard Groovy syntax and the conventions established in `main.nf` and the project's modules.
* **Containerization:** All new processes must be containerized. Ensure your changes reference a public Docker/Singularity image (ideally from StaPHB or a reputable source).

### Submitting a Pull Request (PR)

1.  Push your changes to your fork: `git push origin YOUR-BRANCH-NAME`
2.  Open a **Pull Request** targeting the **`main`** branch of the `UPHL-BioNGS/Cecret` repository.
3.  In the PR description, provide:
    * A **clear title** summarizing the change.
    * A **detailed description** of what you changed and why.
    * A reference to the original issue it addresses (e.g., "Closes #123").
    * A note confirming that you have tested the change.

---

## 🏛️ Code of Conduct

Please note that this project is governed by a **Code of Conduct**. By participating, you are expected to uphold this code. Please report unacceptable behavior to (eriny@utah.gov)[eriny@utah.gov].
