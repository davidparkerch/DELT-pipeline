"""Supervised baseline: cross-validated random forest with scikit-learn."""

from typing import Any

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import KFold, StratifiedKFold, cross_validate

from .._numbers import significant
from ..io import FeatureTable

N_ESTIMATORS = 200
MIN_FOLDS = 2
BINARY = 2


def evaluate_model(table: FeatureTable, cv_folds: int, seed: int, threads: int = 1) -> dict[str, Any] | None:
    """Cross-validate a random forest baseline on the target and rank feature importance.

    The baseline tells downstream ML work how much signal the features carry.

    Returns:
        Metrics (mean and std across folds) and importances, a ``skipped`` reason when the
        data cannot support cross-validation, or ``None`` without a target.
    """
    y, classes = table.encoded_target()
    if y is None:
        return None
    X = table.matrix()  # noqa: N806 - conventional ML name
    if table.task == "classification":
        min_class = int(np.bincount(y).min()) if classes and len(classes) > 1 else 0
        folds = min(cv_folds, min_class)
        estimator: Any = RandomForestClassifier(n_estimators=N_ESTIMATORS, random_state=seed, n_jobs=threads)
        splitter: Any = StratifiedKFold(n_splits=max(folds, MIN_FOLDS), shuffle=True, random_state=seed)
        scoring = {"accuracy": "accuracy", "balanced_accuracy": "balanced_accuracy", "f1_macro": "f1_macro"}
        if classes is not None and len(classes) == BINARY:
            scoring["roc_auc"] = "roc_auc"
    else:
        folds = min(cv_folds, len(y))
        estimator = RandomForestRegressor(n_estimators=N_ESTIMATORS, random_state=seed, n_jobs=threads)
        splitter = KFold(n_splits=max(folds, MIN_FOLDS), shuffle=True, random_state=seed)
        scoring = {"r2": "r2", "mae": "neg_mean_absolute_error", "rmse": "neg_root_mean_squared_error"}

    result: dict[str, Any] = {"task": table.task, "estimator": type(estimator).__name__, "folds": folds}
    if folds < MIN_FOLDS:
        result["skipped"] = "not enough samples per class/fold for cross-validation"
        return result

    scores = cross_validate(estimator, X, y, cv=splitter, scoring=scoring)
    metrics = {}
    for name, scorer in scoring.items():
        values = scores[f"test_{name}"]
        values = -values if scorer.startswith("neg_") else values
        metrics[name] = {"mean": significant(values.mean()), "std": significant(values.std())}
    estimator.fit(X, y)
    ranked = sorted(zip(table.features, estimator.feature_importances_, strict=True), key=lambda item: -item[1])
    result["metrics"] = metrics
    result["feature_importance"] = [
        {"feature": name, "importance": significant(importance)} for name, importance in ranked
    ]
    return result
