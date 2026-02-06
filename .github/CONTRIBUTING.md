* This is a work in progress and not ready for public use, but you are welcome to try it out or help improving it.

* **Code of conduct:** Priority is given to a clean commit history and saving time over user convenience:
    * Commits are pushed directly to `master` and `main` branches.
    * Developers may force-push/rebase recently published commits. 
        * **You may need the CLI to sync forks.**
        * **The commit history of contributions (PRs) will be preserved.**
    * Tags/releases will be added when the code stabilizes.

* LLM-generated code/docs are allowed but require human review.

* Documentation is sparse. Refer to [tests/](../tests/) and code snippets. Help is welcome.

* **Branches:**
    * `main`: Documentation and releases. Expected to pass CI, though currently unstable and not fully tested.
    * `master`: CI/CD testing ([workflows](.github/workflows/)). May be unstable or out-of-sync with `main`.
