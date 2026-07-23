from quantas.core.events import Event, EventLevel, ListObserver


def test_event_progress_is_valid():
    event = Event("Half done", level=EventLevel.PROGRESS, progress=0.5)
    assert event.progress == 0.5


def test_list_observer_collects_events():
    observer = ListObserver()
    observer(Event("Test message"))

    assert len(observer.events) == 1
    assert observer.events[0].message == "Test message"
