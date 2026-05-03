import { useCallback, useEffect, useMemo, useRef, useState } from "react";
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
import { getRandomString } from "../../utilities/Utilities";
import { useChartManager } from "@/react/hooks";
import type { DGEChartManager } from "../../datastore/dgeIntegration";
import { DEFAULT_DGE_BATCH_SIZE } from "../../lib/constants";

interface DGEDialogContentProps {
	dataStore: DataStore;
	onClose: () => void;
}

function DGEDialogContent({ dataStore, onClose }: DGEDialogContentProps) {
	const chartManager = useChartManager() as unknown as DGEChartManager;
	const [groupColumn, setGroupColumn] = useState("");
	const [targetGroup, setTargetGroup] = useState("");
	const [referenceGroup, setReferenceGroup] = useState("rest");
	const [selectedGenesDsName, setSelectedGenesDsName] = useState("");
	const [running, setRunning] = useState(false);
	const [progress, setProgress] = useState({ done: 0, total: 0 });
	const [result, setResult] = useState<DGEIntegrationResult | null>(null);
	const [error, setError] = useState<string | null>(null);
	const isMountedRef = useRef(true);
	const runAbortRef = useRef<AbortController | null>(null);

	useEffect(() => {
		// Keep this resilient in React StrictMode (dev), where mount effects
		// are intentionally mounted/cleaned up/re-mounted.
		isMountedRef.current = true;
		return () => {
			isMountedRef.current = false;
			runAbortRef.current?.abort();
			runAbortRef.current = null;
		};
	}, []);

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

	const safeTargetGroup = selectedColumnValues.includes(targetGroup) ? targetGroup : "";
	const safeReferenceGroup = referenceGroup === "rest" || selectedColumnValues.includes(referenceGroup)
		? referenceGroup : "rest";

	const genesDatasourceOptions = useMemo(() => {
		return findDGECapableDatasources(chartManager)
			.filter((p) => p.cellsDsName === dataStore.name)
			.map((p) => p.genesDsName)
			.sort((a, b) => a.localeCompare(b));
	}, [chartManager, dataStore.name]);

	const referenceOptions = useMemo(() => {
		return ["rest", ...selectedColumnValues.filter((v: string) => v !== safeTargetGroup)];
	}, [selectedColumnValues, safeTargetGroup]);

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

	useEffect(() => {
		if (genesDatasourceOptions.length === 0) {
			setSelectedGenesDsName("");
			return;
		}
		if (!genesDatasourceOptions.includes(selectedGenesDsName)) {
			setSelectedGenesDsName(genesDatasourceOptions[0]);
		}
	}, [genesDatasourceOptions, selectedGenesDsName]);

	const handleAbort = useCallback(() => {
		const currentRun = runAbortRef.current;
		if (!currentRun) return;
		currentRun.abort();
		runAbortRef.current = null;
		if (!isMountedRef.current) return;
		setRunning(false);
		setProgress({ done: 0, total: 0 });
		setError("DGE run aborted.");
	}, []);

	const handleRun = useCallback(async () => {
		const w: any = window as any;
		const getBrowserMemoryBudgetBytes = () => {
			const jsHeapLimit = w?.performance?.memory?.jsHeapSizeLimit;
			const deviceMemGB = w?.navigator?.deviceMemory;
			const candidates: number[] = [];
			if (typeof jsHeapLimit === "number" && Number.isFinite(jsHeapLimit) && jsHeapLimit > 0) {
				candidates.push(jsHeapLimit * 0.6);
			}
			if (typeof deviceMemGB === "number" && Number.isFinite(deviceMemGB) && deviceMemGB > 0) {
				candidates.push(deviceMemGB * 1024 * 1024 * 1024 * 0.25);
			}
			if (!candidates.length) return null;
			return Math.min(...candidates);
		};

		const cellsCount = dataStore.size ?? 0;
		const genesCount = selectedGenesDsName
			? chartManager?.dsIndex?.[selectedGenesDsName]?.dataStore?.size ?? null
			: null;
		const budgetBytes = getBrowserMemoryBudgetBytes();
		// Better estimate: DGE loads gene expression columns in batches.
		// Each gene expression column is a Float32Array over all cells (~4 bytes per cell).
		const batchSize = DEFAULT_DGE_BATCH_SIZE;
		const bytesPerGeneColumn = cellsCount * 4;
		const bytesPerBatch = bytesPerGeneColumn * batchSize;
		const probeColumns = 50; // from runDGEOnDataStore MAX_PROBE
		const bytesForProbe = bytesPerGeneColumn * Math.min(probeColumns, genesCount ?? probeColumns);
		// Shared buffers always present:
		const bytesForFilter = cellsCount; // Uint8 filterBuffer
		const bytesForGroup = cellsCount; // Uint8 groupSAB
		// DGE result columns written to genes DS (6 Float32 arrays sized by gene count).
		const bytesForResults = (genesCount ?? 0) * 4 * 6;
		// Peak-ish estimate: one batch of gene columns + existing buffers + results + some overhead factor.
		const estimatedMinBytes = bytesPerBatch + bytesForProbe + bytesForFilter + bytesForGroup + bytesForResults;
		const headroomFactor = 2.0;
		const budgetUnknown = budgetBytes == null;
		const shouldWarn =
			// If we have a budget, compare estimate to it.
			(!budgetUnknown && estimatedMinBytes * headroomFactor > (budgetBytes as number)) ||
			// If we *don't* have a budget signal, only warn for clearly-large runs.
			(budgetUnknown && (cellsCount >= 500_000 || (genesCount ?? 0) >= 50_000));

		if (shouldWarn) {
			const estimatedMB = estimatedMinBytes / (1024 * 1024);
			const budgetMB = budgetBytes == null ? null : budgetBytes / (1024 * 1024);
			const msg =
				`This DGE run will process ${cellsCount.toLocaleString()} cells` +
				(genesCount ? ` and ${genesCount.toLocaleString()} genes` : "") +
				". This is likely to consume a lot of browser memory and may crash (OOM / Array buffer allocation failed)." +
				`\n\nEstimate (very rough): ~${estimatedMB.toFixed(0)} MB minimum` +
				(budgetMB == null ? "" : `; budget ~${budgetMB.toFixed(0)} MB`) +
				".\n\nContinue?";

			if (!window.confirm(msg)) {
				return;
			}
		}

		runAbortRef.current?.abort();
		const abortController = new AbortController();
		runAbortRef.current = abortController;

		setRunning(true);
		setError(null);
		setResult(null);
		setProgress({ done: 0, total: 0 });

		console.log("=== DGE RUN START ===");
		console.log("groupColumn:", groupColumn, "target:", safeTargetGroup, "reference:", safeReferenceGroup);
		console.log("dataStore:", dataStore.name, "size:", dataStore.size);

		try {
			const cm = chartManager;
			if (!selectedGenesDsName) {
				throw new Error("Select a linked genes datasource");
			}

			const dgeResult = await runDGEOnDataStore(
				cm,
				dataStore.name,
				selectedGenesDsName,
				{
					groupColumn,
					targetGroup: safeTargetGroup,
					referenceGroup: safeReferenceGroup,
				},
				(done, total) => {
					if (isMountedRef.current && !abortController.signal.aborted) {
						setProgress({ done, total });
					}
				},
				{ signal: abortController.signal },
			);
			if (!isMountedRef.current || abortController.signal.aborted) return;

			setResult(dgeResult);

			// Verify DGE columns exist on genes DS before charting
			const gDS = cm.dsIndex[selectedGenesDsName]?.dataStore;
			if (gDS) {
				for (const f of ["dge_effect_size", "dge_neg_log10_pval_adj"]) {
					const c = gDS.columnIndex[f];
					if (!c) { console.error("[DGE] column missing:", f); continue; }
					const d = c.data;
					if (!d) { console.error("[DGE] column has no data:", f); continue; }
					let valid = 0, nanCount = 0, min = Infinity, max = -Infinity;
					for (let i = 0; i < d.length; i++) {
						if (Number.isNaN(d[i])) { nanCount++; } else { valid++; if (d[i] < min) min = d[i]; if (d[i] > max) max = d[i]; }
					}
					console.log(`[DGE] genes DS col "${f}": length=${d.length}, valid=${valid}, NaN=${nanCount}, min=${min}, max=${max}, minMax=`, c.minMax);
				}
				console.log("[DGE] genes DS size:", gDS.size, "columnsWithData:", gDS.columnsWithData.filter((x: string) => x.startsWith("dge_")));
			}

			const volcanoChartId = `dge_volcano_${getRandomString()}`;
			try {
				await cm.addChart(selectedGenesDsName, {
					id: volcanoChartId,
					type: "wgl_scatter_plot",
					title: `DGE: ${safeTargetGroup} vs ${safeReferenceGroup}`,
					param: ["dge_effect_size", "dge_neg_log10_pval_adj"],
					size: [500, 400],
				}, true);
			} catch (err: any) {
				if (abortController.signal.aborted) return;
				console.warn("DGE completed but opening volcano plot failed", {
					error: err,
					volcanoChartId,
					selectedGenesDsName,
					safeTargetGroup,
					safeReferenceGroup,
				});
			}
		} catch (e: any) {
			if (abortController.signal.aborted || !isMountedRef.current) return;
			if (e?.name !== "AbortError") {
				setError(e.message || String(e));
			}
		} finally {
			const isCurrentRun = runAbortRef.current === abortController;
			if (isCurrentRun) {
				runAbortRef.current = null;
			}
			// Always settle "running" for the currently tracked run, including aborts.
			if (isMountedRef.current && isCurrentRun) {
				setRunning(false);
			}
		}
	}, [chartManager, dataStore.name, dataStore.size, groupColumn, safeTargetGroup, safeReferenceGroup, selectedGenesDsName]);

	const sigCount = result
		? result.results.filter((r) => r.pvalAdj < 0.05).length
		: 0;

	return (
		<Box sx={{ p: 2, display: "flex", flexDirection: "column", gap: 2, minWidth: 350 }}>
			<FormControl fullWidth size="small">
				<InputLabel>DGE Datasource</InputLabel>
				<Select
					value={selectedGenesDsName}
					label="Genes Datasource"
					onChange={(e) => setSelectedGenesDsName(e.target.value)}
					disabled={running || genesDatasourceOptions.length <= 1}
				>
					{genesDatasourceOptions.map((name: string) => (
						<MenuItem key={name} value={name}>
							{name}
						</MenuItem>
					))}
				</Select>
			</FormControl>

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
					value={safeTargetGroup}
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
					value={safeReferenceGroup}
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
				disabled={running || !groupColumn || !safeTargetGroup || !selectedGenesDsName}
				fullWidth
			>
				{running ? "Running..." : "Run DGE"}
			</Button>

			{running && (
				<Box sx={{ display: "flex", flexDirection: "column", gap: 1 }}>
					<LinearProgress
						variant={progress.total > 0 ? "determinate" : "indeterminate"}
						value={progress.total > 0 ? (progress.done / progress.total) * 100 : 0}
					/>
					<Typography variant="caption" color="text.secondary">
						{progress.total > 0
							? `Loading gene columns... (batch ${progress.done}/${progress.total})`
							: "Initializing..."}
					</Typography>
					<Button
						variant="outlined"
						color="warning"
						onClick={handleAbort}
						fullWidth
					>
						Abort DGE
					</Button>
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
					{result.skippedGenes > 0 && (
						<Typography variant="caption" display="block" sx={{ mt: 0.5, opacity: 0.85 }}>
							{result.skippedGenes} genes skipped (failed to load)
						</Typography>
					)}
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
