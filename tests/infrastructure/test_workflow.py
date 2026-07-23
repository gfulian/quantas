from quantas.core.events import EventLevel, ListObserver
from quantas.models import BasicCalculator, InputData, ResultData, ResultMetadata


class DummyCalculator(BasicCalculator):
    """
    Small calculator used to test the common Quantas workflow.
    """

    module = "dummy"
    method = "test"

    def run(self) -> ResultData:
        self.emit("Running dummy calculation", level=EventLevel.INFO, progress=0.5)

        value = self.input_data.data["x"] * self.options.get("factor", 1.0)

        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=self.options,
            results={"value": value},
            warnings=self.warnings,
        )


def test_dummy_calculator_workflow():
    observer = ListObserver()
    input_data = InputData(data={"x": 2.0})

    calculator = DummyCalculator(
        input_data=input_data,
        options={"factor": 3.0},
        observer=observer,
    )

    result = calculator.execute()

    assert calculator.completed is True
    assert calculator.error is None
    assert result.results["value"] == 6.0
    assert result.metadata.module == "dummy"

    messages = [event.message for event in observer.events]

    assert "Preparing calculation" in messages
    assert "Running dummy calculation" in messages
    assert "Finalizing calculation" in messages
    assert "Calculation completed" in messages


class ProgressCalculator(BasicCalculator):
    """Calculator used to verify ephemeral progress-event persistence."""

    module = "progress"
    method = "test"

    def run(self) -> ResultData:
        self.emit(
            "Working",
            level=EventLevel.PROGRESS,
            progress=0.5,
            data={"current": 1, "total": 2},
        )
        self.emit("A persistent warning", level=EventLevel.WARNING)
        return ResultData(
            metadata=ResultMetadata(module=self.module, method=self.method),
            input_data=self.input_data,
            options=self.options,
        )


def test_progress_events_are_observed_but_not_persisted() -> None:
    """High-frequency progress updates do not inflate HDF5 workflow logs."""
    observer = ListObserver()
    result = ProgressCalculator(InputData(), observer=observer).execute()
    assert any(event.level is EventLevel.PROGRESS for event in observer.events)
    assert all(record.level != EventLevel.PROGRESS.value for record in result.events)
    assert any(record.level == EventLevel.WARNING.value for record in result.events)
