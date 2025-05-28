#!/bin/bash
# Screen-optimized RNA-seq pipeline runner
# Ensures proper logging and screen session management

set -euo pipefail

# Configuration
SCREEN_NAME="rnaseq_pipeline"
LOG_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SCREEN_LOG="$LOG_DIR/screen_session_${TIMESTAMP}.log"

# Ensure log directory exists
mkdir -p "$LOG_DIR"

# Function to run pipeline in screen
run_in_screen() {
    local force_flag=""
    if [[ "${1:-}" == "--force" ]]; then
        force_flag="--force"
    fi
    
    echo "Starting RNA-seq pipeline in screen session: $SCREEN_NAME"
    echo "Screen log: $SCREEN_LOG"
    echo "Pipeline arguments: $force_flag"
    
    # Kill existing screen session if it exists
    if screen -list | grep -q "$SCREEN_NAME"; then
        echo "Terminating existing screen session..."
        screen -S "$SCREEN_NAME" -X quit 2>/dev/null || true
        sleep 2
    fi
    
    # Start new screen session with pipeline
    screen -dmS "$SCREEN_NAME" -L -Logfile "$SCREEN_LOG" bash -c "
        cd /mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/scripts/pipeline_v2.2
        echo 'Starting pipeline at \$(date)'
        echo 'Working directory: \$(pwd)'
        echo 'Arguments: $force_flag'
        
        # Set environment for better output handling
        export PYTHONUNBUFFERED=1
        export TERM=xterm
        
        # Run the pipeline
        bash run_rnaseq_pipeline.sh $force_flag
        
        echo 'Pipeline completed at \$(date)'
        echo 'Press any key to exit screen session...'
        read -n 1
    "
    
    echo "Pipeline started in screen session: $SCREEN_NAME"
    echo ""
    echo "To monitor the pipeline:"
    echo "  screen -r $SCREEN_NAME"
    echo ""
    echo "To detach from screen (while keeping pipeline running):"
    echo "  Ctrl+A, then D"
    echo ""
    echo "To view the log file:"
    echo "  tail -f $SCREEN_LOG"
    echo ""
    echo "To check if screen session is running:"
    echo "  screen -list"
}

# Function to monitor existing session
monitor_session() {
    if screen -list | grep -q "$SCREEN_NAME"; then
        echo "Attaching to existing screen session: $SCREEN_NAME"
        screen -r "$SCREEN_NAME"
    else
        echo "No active screen session found with name: $SCREEN_NAME"
        echo "Available screen sessions:"
        screen -list
    fi
}

# Function to show session status
show_status() {
    echo "=== Screen Session Status ==="
    if screen -list | grep -q "$SCREEN_NAME"; then
        echo "✅ Pipeline screen session is running: $SCREEN_NAME"
        screen -list | grep "$SCREEN_NAME"
    else
        echo "❌ No pipeline screen session found"
    fi
    
    echo ""
    echo "=== Recent Log Files ==="
    if ls -la "$LOG_DIR"/screen_session_*.log 2>/dev/null; then
        echo "Latest screen log:"
        ls -1t "$LOG_DIR"/screen_session_*.log | head -1
    fi
    
    if ls -la "$LOG_DIR"/pipeline_runs/pipeline_run_*.log 2>/dev/null; then
        echo "Latest pipeline log:"
        ls -1t "$LOG_DIR"/pipeline_runs/pipeline_run_*.log | head -1
    fi
}

# Function to clean up old sessions
cleanup_sessions() {
    echo "Cleaning up old screen sessions..."
    screen -wipe 2>/dev/null || true
    
    echo "Removing old log files (older than 7 days)..."
    find "$LOG_DIR" -name "screen_session_*.log" -mtime +7 -delete 2>/dev/null || true
    find "$LOG_DIR/pipeline_runs" -name "pipeline_run_*.log" -mtime +7 -delete 2>/dev/null || true
}

# Main script logic
case "${1:-}" in
    --monitor|-m)
        monitor_session
        ;;
    --status|-s)
        show_status
        ;;
    --cleanup|-c)
        cleanup_sessions
        ;;
    --force)
        run_in_screen "--force"
        ;;
    --help|-h)
        echo "Usage: $0 [OPTIONS]"
        echo ""
        echo "Options:"
        echo "  (no args)         Start pipeline in screen session"
        echo "  --force           Start pipeline with --force flag"
        echo "  --monitor, -m     Attach to existing screen session"
        echo "  --status, -s      Show status of screen sessions and logs"
        echo "  --cleanup, -c     Clean up old sessions and logs"
        echo "  --help, -h        Show this help message"
        echo ""
        echo "Screen session management:"
        echo "  Screen name: $SCREEN_NAME"
        echo "  Log directory: $LOG_DIR"
        ;;
    *)
        run_in_screen
        ;;
esac