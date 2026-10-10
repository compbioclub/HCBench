__all__ = ["RealBench", "realbench"]


def __getattr__(name):
    if name == "RealBench":
        from .realbench import RealBench

        globals()[name] = RealBench
        return RealBench
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
