import os
import json
import time
import threading
import multiprocessing
from mdvtools.mdvproject import save_json_atomic

def test_json_saving(test_dir="/app/mdv/test_json", iterations=50):
    """
    Test both old and new JSON saving methods under stress.
    Returns comparative results to assess reliability.
    """
    # Recreate old method for comparison
    def save_json_old(file, data):
        with open(file, 'w') as o:
            o.write(json.dumps(data, indent=2, allow_nan=False))
        
    
    os.makedirs(test_dir, exist_ok=True)
    old_file = os.path.join(test_dir, "old_method.json")
    new_file = os.path.join(test_dir, "atomic_method.json")
    
    # Create test data that resembles your application data
    test_data = {
        "all_views": ["view1", "view2", "view3"],
        "state": {
            "permission": "edit",
            "dataloading": {"split": 5, "threads": 2}
        },
        "datasources": [
            {"name": "source1", "columns": [{"field": "col1", "name": "Column 1", "datatype": "text"}] * 50},
            {"name": "source2", "columns": [{"field": "col2", "name": "Column 2", "datatype": "double"}] * 50}
        ] * 10  # Duplicate to make it larger
    }
    
    # Test function for interrupted writes
    def test_interrupted_write(save_func, filename):
        def worker():
            for i in range(iterations):
                # Modify data slightly to ensure it's different each time
                modified_data = test_data.copy()
                modified_data["iteration"] = i
                modified_data["timestamp"] = time.time()
                
                try:
                    save_func(filename, modified_data)
                    time.sleep(0.01)  # Small delay
                except Exception:
                    pass
        
        # Start worker process
        proc = multiprocessing.Process(target=worker)
        proc.start()
        
        # Let it run for a bit
        time.sleep(2)
        
        # Force terminate to simulate crash
        proc.terminate()
        proc.join(1)
        
        # Check if JSON is valid
        try:
            with open(filename, 'r') as f:
                json.load(f)
            return True
        except (json.JSONDecodeError, FileNotFoundError):
            return False
    
    # Test function for concurrent writes
    def test_concurrent_write(save_func, filename, threads=5):
        success_count = 0
        
        for _ in range(iterations // 5):  # Fewer iterations since we have multiple threads
            def thread_worker(thread_id):
                modified_data = test_data.copy()
                modified_data["thread_id"] = thread_id
                modified_data["timestamp"] = time.time()
                
                try:
                    save_func(filename, modified_data)
                except Exception:
                    pass
            
            # Start multiple threads writing to the same file
            threads_list = []
            for i in range(threads):
                t = threading.Thread(target=thread_worker, args=(i,))
                threads_list.append(t)
                t.start()
            
            # Wait for all threads to finish
            for t in threads_list:
                t.join()
            
            # Check if file is valid
            try:
                with open(filename, 'r') as f:
                    json.load(f)
                success_count += 1
            except (json.JSONDecodeError, FileNotFoundError):
                pass
        
        return success_count / (iterations // 5)
    
    # Run tests
    results = {
        "interrupted_write": {
            "old_method": test_interrupted_write(save_json_old, old_file),
            "atomic_method": test_interrupted_write(save_json_atomic, new_file)
        },
        "concurrent_write": {
            "old_method": test_concurrent_write(save_json_old, old_file),
            "atomic_method": test_concurrent_write(save_json_atomic, new_file)
        }
    }
    
    return results

# Run the test
if __name__ == "__main__":
    results = test_json_saving()
    print(json.dumps(results, indent=2))