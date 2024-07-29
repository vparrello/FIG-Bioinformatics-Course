# Installing the BV-BRC Command-Line Interface

During some of the later exercises in this course, you will need to fetch genomic and metagenomic data from [BV-BRC, the "Bacteial and Viral Bioinformatics Center"](https://www.bv-brc.org/) via its "Command-Line Interface" application.
The following instructions will walk you through BV-BRC registration, downloading the CLI app, and logging in to BV-BRC.

1. In order to access the data at BV-BRC, you will first need to register. Please go to
https://www.bv-brc.org/docs/quick_references/registration.html,
and follow the istructions.

2. Once you have registered, you will need to download and install the Command-Line Interface toolkit as documented here:<br>
https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#cli-installation<br>
BV-BRC supports three operating-systems; detailed installation-instructions for each OS may be found here:
    * [macOS](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-macos)
    * [LINUX](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-debian-ubuntu-mint-linux)
    * [Windows](https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#installation-on-windows)

3. Once installed, launch the `BV-BRC.app` under `macOS` and `Windows` by double-clicking on the app's icon. (Under LINUX, the app commands will be installed on your default `PATH`, so just open a new terminal window and start typing commands.)

4. Login to BV-BRC using this command:
```
    p3-login your_BVBRC_username
```
which will prompt you for your BV-BRC password;
once logged in, you can issue any P3-command.
