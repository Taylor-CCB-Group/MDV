import json
import random
import string


class BasePlot:
    def __init__(self, title, plot_type, params, size, position, id=None, **kwargs):
        self.plot_data = {
            "title": title,
            "type": plot_type,
            "param": params,
            "size": size,
            "position": position,
            "id": id if id else self.generate_id(),
        }
        # arbitary key-value pairs can be added to the plot_data
        for key, value in kwargs.items():
            self.plot_data[key] = value

    def generate_id(self):
        """Generate a unique ID for the plot"""
        return str("".join(random.choices(string.ascii_letters, k=6)))

    def set_legend(self, legend):
        self.plot_data["legend"] = legend

    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties

    def set_position(self, position):
        self.position = position

    def to_json(self):
        return json.dumps(self.plot_data, indent=2)
