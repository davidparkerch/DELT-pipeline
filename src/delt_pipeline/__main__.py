"""Allow ``python -m <package>``."""

from .cli import main

raise SystemExit(main())
