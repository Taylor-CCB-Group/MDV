import { types, type IAnyType } from "mobx-state-tree";
import { v4 as uuidv4 } from "uuid";
import Column from "./Column";
/**
 * Helper to generate MST refinement types for min/max
 */
function numberWithConstraints(
  schema: any,
  base: IAnyType
): IAnyType {
  if (schema.minimum !== undefined || schema.maximum !== undefined) {
    const min = schema.minimum !== undefined ? schema.minimum : -Infinity;
    const max = schema.maximum !== undefined ? schema.maximum : Infinity;
    return types.refinement(
      schema.title || "NumberWithConstraints",
      base,
      (value: number) => value >= min && value <= max
    );
  }
  return base;
}

/**
 * Main generator: generate MST type from JSON Schema property
 */
function propertyToMST(
  propSchema: any,
  propName?: string
): IAnyType {
  
  if (
    propSchema.type === "object" && propSchema.title === "__column__"
  ) {
    return Column;
  }
  //ENUM
  if (propSchema.enum) {
    const defaultValue = propSchema.default !== undefined ? propSchema.default : propSchema.enum[0];
    return types.optional(
      types.enumeration(propSchema.title || propName || "Enum", propSchema.enum),
      defaultValue
    );
  }

  // TYPE HANDLING
  switch (propSchema.type) {
    case "string":
      return propSchema.default !== undefined
        ? types.optional(types.string, propSchema.default)
        : types.string;

    case "boolean":
      return propSchema.default !== undefined
        ? types.optional(types.boolean, propSchema.default)
        : types.boolean;

    case "number": {
      const numberType = numberWithConstraints(propSchema, types.number);
      return propSchema.default !== undefined
        ? types.optional(numberType, propSchema.default)
        : numberType;
    }

    case "integer": {
      const intType = numberWithConstraints(propSchema, types.integer);
      return propSchema.default !== undefined
        ? types.optional(intType, propSchema.default)
        : intType;
    }

   case "array": {
        // Tuple-style array: items is an array of schemas
        if (Array.isArray(propSchema.items)) {
            const itemTypes = propSchema.items.map((item: any, idx: any) =>
                propertyToMST(item, propName ? `${propName}Item${idx}` : undefined)
            );
            if (itemTypes.length > 0) {
                if (itemTypes.every((t: any, _: any, arr: any[]) => t === arr[0])) {
                    return types.array(itemTypes[0]);
                } else {
                    return types.array(types.union(...itemTypes));
                }
            }
            // fallback
            return types.array(types.frozen());
        }
        // Regular array: items is a single schema
        const itemType = propertyToMST(propSchema.items, propName ? `${propName}Item` : undefined);
        if (propSchema.default !== undefined) {
            return types.optional(types.array(itemType), propSchema.default);
        }
        return types.array(itemType);
    }

    case "object": {
      // Recursively build nested model
      const props = propSchema.properties || {};
      const required: string[] = propSchema.required || [];
      const mstProps: Record<string, IAnyType> = {};

      for (const [key, childSchema] of Object.entries(props)) {
        const childType = propertyToMST(childSchema as any, key);

        if (required.includes(key)) {
          // Required property: must be present
          mstProps[key] = childType;
        } else if ((childSchema as any).default !== undefined) {
          // Optional with default
          mstProps[key] = types.optional(childType, (childSchema as any).default);
        } else {
          // Optional without default
          mstProps[key] = types.maybe(childType);
        }
      }

      return types.model(propSchema.title || propName || "NestedModel", mstProps);
    }

    default:
      throw new Error(`Unsupported or missing type: ${propSchema.type}`);
  }
}

/**
 * Generate a MobX State Tree model from a JSON schema
 */
export default function jsonSchemaToMST(
  schema: any
): IAnyType {
  if (!schema || schema.type !== "object" || !schema.properties) {
    throw new Error("Schema must be a JSON object with properties.");
  }
  const baseModel = propertyToMST(
    {
      ...schema,
      title:schema.title
    },
    schema.title
  );
  // Only models have .actions, so cast to any to avoid TS error
  let model = (baseModel as any).actions
    ? (baseModel as any).actions((self: any) => ({
        set(path: string, value: any) {
          const keys = path.split(".");
          let target = self;
          for (let i = 0; i < keys.length - 1; i++) {
            target = target[keys[i]];
          }
          target[keys[keys.length - 1]] = value;
        }
      }))
    : baseModel;

  // Add preProcessSnapshot to inject a random id if missing
  model = model.preProcessSnapshot((snapshot: any) => {
    if (!snapshot.id) {
      return { ...snapshot, id: uuidv4() };
    }
    return snapshot;
  });

  return model;
}

/**
 * Example usage:
 * 
 * import schema from "./yourSchema.json";
 * import { jsonSchemaToMST } from "./jsonSchemaToMST";
 * 
 * const ConfigModel = jsonSchemaToMST(schema, "ConfigModel");
 * 
 * // Create instance
 * const config = ConfigModel.create({
 *   theme: "light",
 *   showSidebar: true,
 *   items: [{ label: "First", type: "secondary" }]
 * });
 * 
 * config.theme = "dark"; // reactive!
 */