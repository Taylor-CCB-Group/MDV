import { Box, Typography, Grid, Paper } from '@mui/material';
import H5JsonViewer from './H5JsonViewer';

interface H5MetadataPreviewProps {
  metadata: {
    uns: Record<string, any>;
    obs: Record<string, any>;
    var: Record<string, any>;
  };
}

const H5MetadataPreview = ({ metadata }: H5MetadataPreviewProps) => {
  return (
    <Box sx={{ width: '100%', px: 2 }}>
      <Paper 
        elevation={4}
        sx={{ 
          overflowX: 'auto',
          border: '1px solid',
          borderColor: 'var(--main_panel_color)',
          position: 'relative'  // Added to establish stacking context
        }}
      >
        <Paper
          elevation={1}
          sx={{
            position: 'sticky',
            zIndex: 1,
            borderColor: 'var(--main_panel_color)',
          }}
        >
          <Typography 
            variant="h6" 
            sx={{ 
              p: 2, 
              textAlign: 'center',
              backgroundColor: 'var(--background_color)'
            }}
          >
            H5 File Metadata
          </Typography>
        </Paper>
        
        <Box sx={{ 
          minWidth: 'fit-content',
          maxWidth: '100%',
          p: 2,
          backgroundColor: 'var(--primary_background)',
          position: 'relative',  // For proper stacking
          zIndex: 0  // Below the header
        }}>
          <Grid 
            container 
            spacing={2}
            sx={{
              flexWrap: 'nowrap',
              width: 'auto'
            }}
          >
            <Grid item sx={{ minWidth: '350px' }}>
              <Paper 
                elevation={0} 
                sx={{ 
                  p: 2,
                  backgroundColor: 'white',
                  height: '100%',
                  overflowX: 'auto',
                  border: '1px solid',
                  borderColor: 'rgb(209, 213, 219)'
                }}
              >
                <H5JsonViewer
                  data={metadata.uns}
                  title="Unstructured Data (uns)"
                  initiallyExpanded={true}
                />
              </Paper>
            </Grid>

            <Grid item sx={{ minWidth: '350px' }}>
              <Paper 
                elevation={0} 
                sx={{ 
                  p: 2,
                  backgroundColor: 'white',
                  height: '100%',
                  overflowX: 'auto',
                  border: '1px solid',
                  borderColor: 'rgb(209, 213, 219)'
                }}
              >
                <H5JsonViewer
                  data={metadata.obs}
                  title="Observation Data (obs)"
                  initiallyExpanded={true}
                />
              </Paper>
            </Grid>

            <Grid item sx={{ minWidth: '350px' }}>
              <Paper 
                elevation={0} 
                sx={{ 
                  p: 2,
                  backgroundColor: 'white',
                  height: '100%',
                  overflowX: 'auto',
                  border: '1px solid',
                  borderColor: 'rgb(209, 213, 219)'
                }}
              >
                <H5JsonViewer
                  data={metadata.var}
                  title="Variable Data (var)"
                  initiallyExpanded={true}
                />
              </Paper>
            </Grid>
          </Grid>
        </Box>
      </Paper>
    </Box>
  );
};

export default H5MetadataPreview;