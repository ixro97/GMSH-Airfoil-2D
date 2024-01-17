class InstanceTracker(type):
    def __init__(cls, name, bases, attrs):
        super().__init__(name, bases, attrs)
        cls.instances = []

    def __call__(cls, *args, **kwargs):
        instance = super().__call__(*args, **kwargs)
        cls.instances.append(instance)
        return instance