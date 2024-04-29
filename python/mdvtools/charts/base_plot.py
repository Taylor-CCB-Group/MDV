import json

class BasePlot:
    def __init__(self, title, plot_type, params, size, position, id=None):
        self.plot_data = {
            "title": title,
            "type": plot_type,
            "param": params,
            "size": size,
            "position": position,
            "id": id if id else self.generate_id()
        }

    def generate_id(self):
        # Generate a unique ID for the plot
        return "unique_id"
        # return str(uuid.uuid4()) ?


    def set_legend(self, legend):
        self.plot_data["legend"] = legend

    def set_axis_properties(self, axis, properties):
        if "axis" not in self.plot_data:
            self.plot_data["axis"] = {}
        self.plot_data["axis"][axis] = properties

    def to_json(self):
        return json.dumps(self.plot_data, indent=2)
