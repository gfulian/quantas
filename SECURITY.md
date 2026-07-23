# Security policy

## Supported versions

Quantas 2 is currently in beta.  Security fixes are applied to the latest published
beta or release candidate.  The legacy Quantas 1 series is maintained only where a
separate public release statement explicitly says so.

## Reporting a vulnerability

Do not open a public issue for a vulnerability that could expose user data, execute
untrusted code, corrupt scientific results, or compromise a release artifact.  Use
GitHub's private vulnerability-reporting facility for the Quantas repository when it
is enabled.  If that facility is unavailable, contact the maintainers through the
address published in the repository metadata and clearly mark the message as a
security report.

Include:

- affected Quantas version and platform;
- minimal reproduction steps;
- expected and observed impact;
- whether untrusted input is required;
- any proposed mitigation.

Scientific disagreements, numerical accuracy reports, malformed external-code
outputs, and ordinary crashes should use the normal issue tracker unless they also
create a security impact.
