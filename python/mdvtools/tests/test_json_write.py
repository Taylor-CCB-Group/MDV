import os
import json
import tempfile
import time
import threading
import random
import multiprocessing
import shutil
from pathlib import Path
from mdvtools.mdvproject import save_json_atomic

def save_json_old(file, data):
    """Recreate the old implementation for testing"""
    with open(file, 'w') as o:
        o.write(json.dumps(data, indent=2, allow_nan=False))
        o.close()

class JSONTruncationTests:
    def __init__(self, test_dir="/app/mdv/json_test"):
        self.test_dir = Path(test_dir)
        self.test_dir.mkdir(exist_ok=True)
        
    def cleanup(self):
        """Clean up test files after running tests"""
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
            
    def generate_large_json(self, size=1000):
        """Generate a large nested JSON object"""
        data = {
            "metadata": {
                "id": "test-large-json",
                "description": "A large JSON object for testing file writes"
            },
            "items": []
        }
        
        # Add many nested items to increase JSON size
        for i in range(size):
            data["items"].append({
                "id": f"item-{i}",
                "properties": {
                    "name": f"Test Item {i}",
                    "value": random.random(),
                    "tags": [f"tag-{j}" for j in range(10)],
                    "nested": {
                        "level1": {
                            "level2": {
                                "level3": f"Deep nested value {i}"
                            }
                        }
                    }
                }
            })
            
        return data
    
    def validate_json_file(self, path):
        """Check if a JSON file is valid and complete"""
        try:
            with open(path, 'r') as f:
                data = json.load(f)
                # Additional validation can be added here
                return True
        except (json.JSONDecodeError, FileNotFoundError):
            return False
        
    def interrupted_write_test(self, use_atomic=False):
        """Test where the process is killed during JSON writing"""
        test_file = self.test_dir / "interrupted.json"
        
        # Create a child process to perform the write
        def child_process():
            large_data = self.generate_large_json(5000)  # Make it large enough to take time
            save_func = save_json_atomic if use_atomic else save_json_old
            
            # Signal that we're starting to write
            with open(self.test_dir / "started.flag", "w") as f:
                f.write("1")
                
            # Perform the slow write operation
            save_func(str(test_file), large_data)
            
            # Signal completion
            with open(self.test_dir / "completed.flag", "w") as f:
                f.write("1")
                
        # Start the child process
        proc = multiprocessing.Process(target=child_process)
        proc.start()
        
        # Wait for child to start writing
        start_time = time.time()
        while not (self.test_dir / "started.flag").exists():
            time.sleep(0.01)
            if time.time() - start_time > 5:  # 5 second timeout
                proc.terminate()
                return "TIMEOUT_WAITING_TO_START"
        
        # Give it a moment to get into the actual writing
        time.sleep(0.1)
        
        # Kill the process while it's writing
        proc.terminate()
        
        # Check if the file exists and is valid
        if not test_file.exists():
            return "FILE_NOT_CREATED"
        
        # Validate the JSON
        if self.validate_json_file(test_file):
            return "VALID_JSON"
        else:
            return "INVALID_JSON"
    
    def concurrent_write_test(self, use_atomic=False, iterations=20):
        """Test multiple processes writing to the same file"""
        test_file = self.test_dir / "concurrent.json"
        success_count = 0
        
        def worker(worker_id):
            data = {
                "worker_id": worker_id,
                "timestamp": time.time(),
                "data": self.generate_large_json(100)
            }
            
            save_func = save_json_atomic if use_atomic else save_json_old
            try:
                save_func(str(test_file), data)
                return True
            except Exception:
                return False
            
        for i in range(iterations):
            # Run multiple workers concurrently
            threads = []
            for j in range(5):  # 5 concurrent writers
                thread = threading.Thread(target=lambda: worker(f"{i}-{j}"))
                threads.append(thread)
                thread.start()
                
            # Wait for all threads to complete
            for thread in threads:
                thread.join()
                
            # Check if the file is valid
            if self.validate_json_file(test_file):
                success_count += 1
                
            # Small delay between iterations
            time.sleep(0.01)
            
        return success_count / iterations
    
    def rapid_succession_test(self, use_atomic=False, iterations=100):
        """Test rapid file writes in succession"""
        test_file = self.test_dir / "rapid.json"
        success_count = 0
        
        for i in range(iterations):
            data = {
                "iteration": i,
                "timestamp": time.time(),
                "random_data": [random.random() for _ in range(100)]
            }
            
            save_func = save_json_atomic if use_atomic else save_json_old
            try:
                save_func(str(test_file), data)
                
                # Verify file integrity
                if self.validate_json_file(test_file):
                    success_count += 1
            except Exception:
                pass
                
            # No delay to maximize stress
            
        return success_count / iterations
    
    def resource_constrained_test(self, use_atomic=False):
        """Test writing while system resources are constrained"""
        test_file = self.test_dir / "resource.json"
        success_count = 0
        
        # Create resource competition with many parallel operations
        def resource_hog():
            # Create CPU and I/O pressure
            temp_files = []
            for _ in range(10):
                tf = tempfile.NamedTemporaryFile(dir=self.test_dir, delete=False)
                temp_files.append(tf.name)
                # Write random data
                for _ in range(1000):
                    tf.write(os.urandom(1024))
                tf.flush()
                tf.close()
            
            # CPU pressure
            _ = [random.random() ** random.random() for _ in range(100000)]
            
            # Clean up
            for tf in temp_files:
                os.unlink(tf)
        
        for i in range(20):
            # Start resource hogs
            hog_threads = [threading.Thread(target=resource_hog) for _ in range(3)]
            for t in hog_threads:
                t.start()
                
            # Perform the save operation
            data = self.generate_large_json(200)
            save_func = save_json_atomic if use_atomic else save_json_old
            
            try:
                save_func(str(test_file), data)
                
                # Verify file integrity
                if self.validate_json_file(test_file):
                    success_count += 1
            except Exception:
                pass
                
            # Wait for hog threads to complete
            for t in hog_threads:
                t.join()
                
        return success_count / 20
        
    def run_all_tests(self):
        """Run all tests and return results"""
        results = {
            "interrupted_write": {
                "old": self.interrupted_write_test(use_atomic=False),
                "atomic": self.interrupted_write_test(use_atomic=True)
            },
            "concurrent_write": {
                "old": self.concurrent_write_test(use_atomic=False),
                "atomic": self.concurrent_write_test(use_atomic=True)
            },
            "rapid_succession": {
                "old": self.rapid_succession_test(use_atomic=False),
                "atomic": self.rapid_succession_test(use_atomic=True)
            },
            "resource_constrained": {
                "old": self.resource_constrained_test(use_atomic=False),
                "atomic": self.resource_constrained_test(use_atomic=True)
            }
        }
        
        self.cleanup()
        return results

# Example usage
def run_tests():
    tester = JSONTruncationTests()
    results = tester.run_all_tests()
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    run_tests()