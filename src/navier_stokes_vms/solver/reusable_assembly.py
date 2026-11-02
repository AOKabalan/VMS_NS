"""Assemblers for forms evaluated many times, preprocessed once instead of on every call."""

from collections.abc import Callable

import ufl
from firedrake.assemble import AbstractFormAssembler, BaseFormAssembler, get_assembler
from ufl.form import BaseForm, Form


def reusable_assembler(form: BaseForm) -> AbstractFormAssembler:
    """The assembler `fd.assemble` builds, after the same preprocessing, for repeated use.

    `fd.assemble` preprocesses the form again on every call, which dominates the cost of small
    problems with many outer iterations and local sweeps.
    """
    preprocessed = BaseFormAssembler.preprocess_base_form(form)
    return get_assembler(preprocessed, is_base_form_preprocessed=True)


def scalar_functional(form: Form) -> Callable[[], float]:
    """Preprocess a repeated scalar functional, including identically vanishing forms."""
    expanded = ufl.algorithms.expand_derivatives(form)
    if not any(integral.integrand() != 0 for integral in expanded.integrals()):
        return lambda: 0.0
    return reusable_assembler(form).assemble
