import { useCallback, useEffect, useMemo, useState } from "react";
import {
	Box,
	Button,
	FormControl,
	InputLabel,
	LinearProgress,
	MenuItem,
	Select,
	Typography,
} from "@mui/material";
import type DataStore from "../../datastore/DataStore.js";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import {
	runDGEOnDataStore,
	findDGECapableDatasources,
} from "../../datastore/dgeIntegration";
import type { DGEIntegrationResult } from "../../datastore/dgeIntegration";

interface DGEDialogContentProps {
	dataStore: DataStore;
	onClose: () => void;
}

function DGEDialogContent({ dataStore, onClose }: DGEDialogContentProps) {
	const [groupColumn, setGroupColumn] = useState("");
	const [targetGroup, setTargetGroup] = useState("");
	const [referenceGroup, setReferenceGroup] = useState("rest");
	const [running, setRunning] = useState(false);
	const [progress, setProgress] = useState({ done: 0, total: 0 });
	const [result, setResult] = useState<DGEIntegrationResult | null>(null);
	const [error, setError] = useState<string | null>(null);

	const categoricalColumns = useMemo(() => {
		return dataStore.columns.filter(
			(c: any) =>
				(c.datatype === "text" || c.datatype === "text16") && c.values,
		);
	}, [dataStore]);

	const selectedColumnValues = useMemo(() => {
		if (!groupColumn) return [];
		const col = dataStore.columnIndex[groupColumn];
		return col?.values ?? [];
	}, [dataStore, groupColumn]);

	const referenceOptions = useMemo(() => {
		return ["rest", ...selectedColumnValues.filter((v: string) => v !== targetGroup)];
	}, [selectedColumnValues, targetGroup]);

	useEffect(() => {
		if (categoricalColumns.length > 0 && !groupColumn) {
			setGroupColumn(categoricalColumns[0].field);
		}
	}, [categoricalColumns, groupColumn]);

	useEffect(() => {
		if (selectedColumnValues.length > 0) {
			setTargetGroup(selectedColumnValues[0]);
			setReferenceGroup("rest");
		}
	}, [selectedColumnValues]);

	const handleRun = useCallback(async () => {
		setRunning(true);
		setError(null);
		setResult(null);
		setProgress({ done: 0, total: 0 });

		try {
			const cm = window.mdv?.chartManager;
			if (!cm) throw new Error("ChartManager not available");

			const pairs = findDGECapableDatasources(cm);
			const pair = pairs.find((p) => p.cellsDsName === dataStore.name);
			if (!pair) throw new Error("No gene expression link found for this datasource");

			const dgeResult = await runDGEOnDataStore(
				cm,
				pair.cellsDsName,
				pair.genesDsName,
				{
					groupColumn,
					targetGroup,
					referenceGroup,
				},
				(done, total) => setProgress({ done, total }),
			);

			setResult(dgeResult);

			try {
				await cm.addChart(pair.genesDsName, {
					type: "wgl_scatter_plot",
					title: `DGE: ${targetGroup} vs ${referenceGroup}`,
					param: ["dge_effect_size", "dge_neg_log10_pval_adj"],
					size: [500, 400],
				}, true);
			} catch {
				// Volcano plot is optional; don't fail the whole operation
			}
		} catch (e: any) {
			setError(e.message || String(e));
		} finally {
			setRunning(false);
		}
	}, [dataStore.name, groupColumn, targetGroup, referenceGroup]);

	const sigCount = result
		? result.results.filter((r) => r.pvalAdj < 0.05).length
		: 0;

	return (
		<Box sx={{ p: 2, display: "flex", flexDirection: "column", gap: 2, minWidth: 350 }}>
			<FormControl fullWidth size="small">
				<InputLabel>Grouping Column</InputLabel>
				<Select
					value={groupColumn}
					label="Grouping Column"
					onChange={(e) => setGroupColumn(e.target.value)}
					disabled={running}
				>
					{categoricalColumns.map((c: any) => (
						<MenuItem key={c.field} value={c.field}>
							{c.name}
						</MenuItem>
					))}
				</Select>
			</FormControl>

			<FormControl fullWidth size="small">
				<InputLabel>Target Group</InputLabel>
				<Select
					value={targetGroup}
					label="Target Group"
					onChange={(e) => setTargetGroup(e.target.value)}
					disabled={running || selectedColumnValues.length === 0}
				>
					{selectedColumnValues.map((v: string) => (
						<MenuItem key={v} value={v}>
							{v}
						</MenuItem>
					))}
				</Select>
			</FormControl>

			<FormControl fullWidth size="small">
				<InputLabel>Reference</InputLabel>
				<Select
					value={referenceGroup}
					label="Reference"
					onChange={(e) => setReferenceGroup(e.target.value)}
					disabled={running || selectedColumnValues.length === 0}
				>
					{referenceOptions.map((v: string) => (
						<MenuItem key={v} value={v}>
							{v === "rest" ? "Rest (all other groups)" : v}
						</MenuItem>
					))}
				</Select>
			</FormControl>

			<Button
				variant="contained"
				onClick={handleRun}
				disabled={running || !groupColumn || !targetGroup}
				fullWidth
			>
				{running ? "Running..." : "Run DGE"}
			</Button>

			{running && (
				<Box>
					<LinearProgress
						variant={progress.total > 0 ? "determinate" : "indeterminate"}
						value={progress.total > 0 ? (progress.done / progress.total) * 100 : 0}
					/>
					<Typography variant="caption" color="text.secondary">
						{progress.total > 0
							? `Loading gene columns... (batch ${progress.done}/${progress.total})`
							: "Initializing..."}
					</Typography>
				</Box>
			)}

			{error && (
				<Typography color="error" variant="body2">
					{error}
				</Typography>
			)}

			{result && (
				<Box sx={{ p: 1, bgcolor: "success.main", color: "success.contrastText", borderRadius: 1 }}>
					<Typography variant="body2">
						Found {sigCount} significant genes (padj &lt; 0.05) in{" "}
						{(result.elapsed / 1000).toFixed(1)}s
					</Typography>
					<Typography variant="caption">
						Effect size: {result.effectSizeLabel === "log2fc" ? "log2 fold change" : "mean difference"}
					</Typography>
				</Box>
			)}
		</Box>
	);
}

export default class DGEDialogReact extends BaseDialog {
	root: ReturnType<typeof createMdvPortal>;
	constructor(dataStore: DataStore) {
		super(
			{
				title: "Differential Gene Expression",
				width: 400,
				height: 480,
			},
			null,
		);
		this.root = createMdvPortal(
			<DGEDialogContent dataStore={dataStore} onClose={() => this.close()} />,
			this.dialog,
			this,
		);
	}
	close(): void {
		super.close();
		this.root.unmount();
	}
}
