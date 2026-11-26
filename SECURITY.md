# Security policy

## Supported versions

Only the latest version of fasttrop gets security fixes. The latest
version is the version on the `main` branch.

## Report a vulnerability

Do not open a public issue for a vulnerability. Use the private report
of GitHub instead:

1.  Open the page
    <https://github.com/farach/fasttrop/security/advisories/new>.
2.  Describe the problem and the steps to cause it.
3.  Give the version of fasttrop and the operating system.

The maintainer reads each report and replies as soon as possible. The
package has C code, so a report of a memory error is in scope. Examples
are a buffer that is too small, or a read outside of an array.

## Keep private data out of public text

Issues, pull requests, and examples are public. Do not put private data
in them. Use a small panel from
[`trop_simulate()`](https://farach.github.io/fasttrop/reference/trop_simulate.md)
or a panel of random numbers to show a problem. Do not paste
credentials, file paths of your computer, or data that you do not have
the right to share.
