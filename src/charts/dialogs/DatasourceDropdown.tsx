import type React from "react";
import { useState, useEffect } from "react";

export interface DropdownProps {
    options: string[];
    onSelect: (value: string) => void;
}

export const DatasourceDropdown: React.FC<DropdownProps> = ({
    options,
    onSelect,
}) => {
    // Initialize selectedOption with the first option in the array (if it exists)
    const [selectedOption, setSelectedOption] = useState<string | null>(null);

    useEffect(() => {
        // When the component mounts or when options change, select the first option
        if (options.length > 0 && selectedOption !== options[0]) {
            setSelectedOption(options[0]);
            onSelect(options[0]);
        }
    }, [options, selectedOption, onSelect]);

    const handleSelect = (event: React.ChangeEvent<HTMLSelectElement>) => {
        const value = event.target.value;
        setSelectedOption(value);
        onSelect(value);
    };

    return (
        <div className="w-full max-w-xs mx-auto">
            <p className="text-lg text-gray-700 dark:text-white my-1">
                <strong>Select Datasource:</strong>
            </p>
            <select
                value={selectedOption || ""}
                onChange={handleSelect}
                className="block w-full px-3 py-2 border border-gray-300 bg-white rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 sm:text-sm"
            >
                <option value="" disabled>
                    -- Select an option --
                </option>
                {options.map((option, index) => (
                    // biome-ignore lint/suspicious/noArrayIndexKey: not a big deal in this case... would be good to fix ideally
                    <option key={index} value={option}>
                        {option}
                    </option>
                ))}
            </select>
        </div>
    );
};
