from quantas.core.events import Event, EventLevel, ListObserver


def test_event_progress_is_valid():
    event = Event("Half done", level=EventLevel.PROGRESS, progress=0.5)
    assert event.progress == 0.5


def test_list_observer_collects_events():
    observer = ListObserver()
    observer(Event("Test message"))

    assert len(observer.events) == 1
    assert observer.events[0].message == "Test message"


def test_list_observer_uses_identity_semantics_and_independent_storage():
    first = ListObserver()
    second = ListObserver()

    first(Event("First"))

    assert first is not second
    assert first != second
    assert [event.message for event in first.events] == ["First"]
    assert second.events == []
