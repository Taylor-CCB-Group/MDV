import { BaseDialog } from "@/utilities/Dialog";
import { createMdvPortal } from "@/react/react_utils";
import type { ManageGateDialogType } from "./ManageGateDialogComponent";
import ManageGateDialogContent from "./ManageGateDialogComponent";
import { ChartProvider } from "../context";
import BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";

export type ManageGateDialogWrapperProps<T extends BaseConfig> = Omit<ManageGateDialogType, "open"> & {
    chart: BaseChart<T>;
};

class ManageGateDialogWrapper<T extends BaseConfig> extends BaseDialog {
    declare root: ReturnType<typeof createMdvPortal>;
    
    constructor(props: ManageGateDialogWrapperProps<T>) {
        const { chart, ...dialogProps } = props;
        const typeName = BaseChart.types[chart.config.type].name;
        const name = `${chart.config.title} (${typeName})`;
        const config = {
            title: `Manage Gates in "${name}"`,
            width: 500,
            onclose: () => dialogProps.onClose(),
        };
        super(config, props);

        this.root = createMdvPortal(
            <ChartProvider chart={chart}>
                <ManageGateDialogContent
                    {...dialogProps}
                    open={true}
                    onClose={() => {
                        dialogProps.onClose();
                        this.close();
                    }}
                />
            </ChartProvider>,
            this.dialog,
            this,
        );
    }

    close(): void {
        this.root?.unmount();
        super.close();
    }
}

export default ManageGateDialogWrapper;
